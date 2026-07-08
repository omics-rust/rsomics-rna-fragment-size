//! Per-transcript mRNA fragment-size distribution for paired-end RNA-seq,
//! matching RSeQC `RNA_fragment_size.py`. For each BED12 transcript it scans the
//! `read1` alignments starting inside the transcript span whose mate also starts
//! inside it, and reports fragment count, mean, median and population standard
//! deviation.
//!
//! ## Origin
//!
//! Independent Rust reimplementation of RSeQC `RNA_fragment_size.py`
//! (LGPL-2.1+), based on the tool's documented behaviour, the SAM/BAM and BED12
//! format specs, and black-box testing against the upstream binary. No upstream
//! source was read; every rule below was fixed by observing RSeQC's output on
//! crafted alignments.
//!
//! Each fragment is measured from the `read1` mate only (flag `0x40`), using the
//! read's start, the mate start (`PNEXT`), and read1's query-aligned length:
//!
//! ```text
//! size = |mrna(mate_start) - mrna(read_start)| + junctions_crossed + qlen(read1)
//! ```
//!
//! where `mrna(g)` is the count of exonic bases before genomic position `g` (a
//! position in an intron collapses to the start of the following exon),
//! `junctions_crossed` is the difference of the two positions' exon indices, and
//! `qlen(read1)` sums the query-consuming CIGAR ops `M/I/=/X` (soft-clips and
//! deletions excluded, insertions included). RSeQC's `+junctions_crossed` term
//! reproduces a one-base-per-junction offset visible in its output. A read1 is
//! used for a transcript when it is uniquely mapped (mapq ≥ `-q`), its mate is
//! mapped (no proper-pair flag required), `tx_start ≤ read_start < tx_end`, and
//! `tx_start ≤ mate_start ≤ tx_end`; the mate's own alignment may extend past the
//! exons or the transcript end. Transcripts below `-n` fragments print
//! `count 0 0 0`; the reported stdev is population (÷n).
//!
//! License: MIT OR Apache-2.0.
//! Upstream credit: RSeQC <https://rseqc.sourceforge.net/> (LGPL-2.1+).

use std::collections::HashMap;
use std::io::Write;
use std::num::NonZero;
use std::path::Path;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};

// Reads that never contribute a fragment: unmapped, secondary, qcfail,
// duplicate, supplementary.
const SKIP: u16 = 0x0004 | 0x0100 | 0x0200 | 0x0400 | 0x0800;
const PAIRED: u16 = 0x0001;
const MATE_UNMAPPED: u16 = 0x0008;
const READ1: u16 = 0x0040;

struct Transcript {
    chrom: String,
    tx_start: i64,
    tx_end: i64,
    symbol: String,
    exons: Vec<(i64, i64)>,
    cum: Vec<i64>,
}

impl Transcript {
    /// mRNA offset (count of exonic bases before genomic position `g`) paired
    /// with the exon index used to count junctions. A position inside an intron
    /// collapses to the start of the following exon; a position past the last
    /// exon collapses to the transcript's full length.
    fn locate(&self, g: i64) -> (i64, usize) {
        for (i, &(s, e)) in self.exons.iter().enumerate() {
            if g < s {
                return (self.cum[i], i);
            }
            if g < e {
                return (self.cum[i] + (g - s), i);
            }
        }
        let last = self.exons.len() - 1;
        let (s, e) = self.exons[last];
        (self.cum[last] + (e - s), last)
    }
}

/// Per chromosome: `(transcripts sorted as (tx_start, index), widest transcript)`.
type ChromIndex = HashMap<String, (Vec<(i64, usize)>, i64)>;

pub struct FragmentSizes {
    transcripts: Vec<Transcript>,
    /// Collected fragment sizes per transcript, in BED order.
    sizes: Vec<Vec<i64>>,
}

pub fn fragment_sizes(
    bam: &Path,
    bed: &Path,
    mapq: u8,
    workers: NonZero<usize>,
) -> Result<FragmentSizes> {
    let transcripts = parse_bed12(bed)?;
    let mut sizes = vec![Vec::new(); transcripts.len()];

    // Per chromosome (upper-cased to match the header keys): transcripts sorted
    // by start, plus the widest transcript so a position query can bound how far
    // left it must scan.
    let mut by_chrom: ChromIndex = HashMap::new();
    for (i, t) in transcripts.iter().enumerate() {
        let entry = by_chrom
            .entry(t.chrom.to_uppercase())
            .or_insert((Vec::new(), 0));
        entry.0.push((t.tx_start, i));
        entry.1 = entry.1.max(t.tx_end - t.tx_start);
    }
    for (v, _) in by_chrom.values_mut() {
        v.sort_unstable();
    }

    let mut reader = rsomics_bamio::open_with_workers(bam, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).to_uppercase())
        .collect();

    let mut rec = RawRecord::default();
    loop {
        if raw::read_record(reader.get_mut(), &mut rec)? == 0 {
            break;
        }
        let flags = rec.flags();
        if flags & SKIP != 0
            || flags & PAIRED == 0
            || flags & READ1 == 0
            || flags & MATE_UNMAPPED != 0
            || rec.mapping_quality() < mapq
        {
            continue;
        }
        let tid = rec.reference_sequence_id();
        if tid < 0 {
            continue;
        }
        let read_start = i64::from(rec.alignment_start());
        let mate_start = i64::from(rec.mate_alignment_start());
        let qlen = query_len(rec.cigar_ops());
        record_fragment(
            &transcripts,
            &by_chrom,
            &ref_names,
            tid,
            read_start,
            mate_start,
            qlen,
            &mut sizes,
        );
    }
    Ok(FragmentSizes { transcripts, sizes })
}

#[allow(clippy::too_many_arguments)]
fn record_fragment(
    transcripts: &[Transcript],
    by_chrom: &ChromIndex,
    ref_names: &[String],
    tid: i32,
    read_start: i64,
    mate_start: i64,
    qlen: i64,
    sizes: &mut [Vec<i64>],
) {
    let Some(chrom) = ref_names.get(tid as usize) else {
        return;
    };
    let Some((entries, max_span)) = by_chrom.get(chrom) else {
        return;
    };
    let hi = entries.partition_point(|&(s, _)| s <= read_start);
    for &(tx_start, idx) in entries[..hi].iter().rev() {
        if tx_start < read_start - *max_span {
            break;
        }
        let t = &transcripts[idx];
        if read_start >= t.tx_end || mate_start < t.tx_start || mate_start > t.tx_end {
            continue;
        }
        let (a_off, a_ex) = t.locate(read_start);
        let (b_off, b_ex) = t.locate(mate_start);
        let junctions = a_ex.abs_diff(b_ex) as i64;
        sizes[idx].push((a_off - b_off).abs() + junctions + qlen);
    }
}

/// Query-aligned length: query-consuming CIGAR ops that are not soft-clips
/// (`M`, `I`, `=`, `X`). Matches pysam `query_alignment_length`.
fn query_len<I: Iterator<Item = (u8, u32)>>(ops: I) -> i64 {
    let mut len = 0;
    for (op, n) in ops {
        if matches!(op, 0 | 1 | 7 | 8) {
            len += i64::from(n);
        }
    }
    len
}

impl FragmentSizes {
    pub fn write_report<W: Write>(&self, mut w: W, min_frags: usize) -> std::io::Result<()> {
        writeln!(
            w,
            "chrom\ttx_start\ttx_end\tsymbol\tfrag_count\tfrag_mean\tfrag_median\tfrag_std"
        )?;
        for (t, v) in self.transcripts.iter().zip(&self.sizes) {
            write!(
                w,
                "{}\t{}\t{}\t{}\t{}\t",
                t.chrom,
                t.tx_start,
                t.tx_end,
                t.symbol,
                v.len()
            )?;
            if v.len() >= min_frags && !v.is_empty() {
                let (mean, median, std) = stats(v);
                writeln!(
                    w,
                    "{}\t{}\t{}",
                    pyfloat(mean),
                    pyfloat(median),
                    pyfloat(std)
                )?;
            } else {
                writeln!(w, "0\t0\t0")?;
            }
        }
        Ok(())
    }
}

fn stats(v: &[i64]) -> (f64, f64, f64) {
    let n = v.len() as f64;
    let mean = v.iter().sum::<i64>() as f64 / n;
    let mut sorted = v.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    let median = if sorted.len() % 2 == 1 {
        sorted[mid] as f64
    } else {
        (sorted[mid - 1] + sorted[mid]) as f64 / 2.0
    };
    let var = v.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / n;
    (mean, median, var.sqrt())
}

fn parse_bed12(path: &Path) -> Result<Vec<Transcript>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| RsomicsError::Io(std::io::Error::other(format!("reading BED12: {e}"))))?;
    let mut out = Vec::new();
    for line in content.lines() {
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 12 {
            eprintln!("[NOTE: input bed must be 12-column] skipped this line: {line}");
            continue;
        }
        let (Ok(tx_start), Ok(tx_end)) = (c[1].parse::<i64>(), c[2].parse::<i64>()) else {
            continue;
        };
        let sizes = parse_list(c[10]);
        let starts = parse_list(c[11]);
        let exons: Vec<(i64, i64)> = starts
            .iter()
            .zip(sizes.iter())
            .map(|(&s, &sz)| (tx_start + s, tx_start + s + sz))
            .collect();
        let mut cum = Vec::with_capacity(exons.len());
        let mut acc = 0;
        for &(s, e) in &exons {
            cum.push(acc);
            acc += e - s;
        }
        out.push(Transcript {
            chrom: c[0].to_string(),
            tx_start,
            tx_end,
            symbol: c[3].to_string(),
            exons,
            cum,
        });
    }
    Ok(out)
}

/// Format like Python's `str(float)`: integral values keep a trailing `.0`,
/// others use the shortest round-tripping representation.
fn pyfloat(x: f64) -> String {
    if x.fract() == 0.0 {
        format!("{x:.1}")
    } else {
        format!("{x}")
    }
}

fn parse_list(field: &str) -> Vec<i64> {
    field
        .trim_end_matches(',')
        .split(',')
        .filter(|s| !s.is_empty())
        .filter_map(|s| s.parse().ok())
        .collect()
}
