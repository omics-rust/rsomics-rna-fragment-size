//! Per-transcript mRNA fragment-size distribution for paired-end RNA-seq,
//! matching RSeQC `RNA_fragment_size.py`. For each BED12 transcript it pairs
//! uniquely-mapped mates whose alignments both fall inside the transcript's
//! exons and measures the fragment length in transcript (mRNA) coordinates —
//! the genomic span with introns removed — then reports count, mean, median and
//! population standard deviation.
//!
//! ## Origin
//!
//! Independent Rust reimplementation of RSeQC `RNA_fragment_size.py`
//! (LGPL-2.1+), based on the tool's documented behaviour, the SAM/BAM and BED12
//! format specs, and black-box testing against the upstream binary. The exact
//! rules — fragment size = `mRNA(rightmost end) - mRNA(leftmost start)` plus one
//! per exon junction the fragment crosses, both mates required uniquely mapped
//! (mapq ≥ `-q`, no proper-pair flag required), both mates inside the transcript
//! exons, population stdev, and transcripts below `-n` printed as `count 0 0 0`
//! — were established by black-box probing. No upstream source was read.
//!
//! License: MIT OR Apache-2.0.
//! Upstream credit: RSeQC <https://rseqc.sourceforge.net/> (LGPL-2.1+).

use std::collections::HashMap;
use std::io::Write;
use std::num::NonZero;
use std::path::Path;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};

// Reads that never contribute a fragment.
const SKIP: u16 = 0x0004 | 0x0100 | 0x0200 | 0x0400 | 0x0800;
const PAIRED: u16 = 0x0001;

struct Transcript {
    chrom: String,
    tx_start: i64,
    tx_end: i64,
    symbol: String,
    exons: Vec<(i64, i64)>,
    cum: Vec<i64>,
}

impl Transcript {
    /// mRNA offset of a genomic coordinate and the index of the exon it falls in,
    /// or `None` if it lies in an intron. Accepts an exclusive end coordinate
    /// sitting on an exon's right edge.
    fn mrna(&self, g: i64) -> Option<(i64, usize)> {
        for (i, &(s, e)) in self.exons.iter().enumerate() {
            if g >= s && g <= e {
                return Some((self.cum[i] + (g - s), i));
            }
        }
        None
    }

    /// True when the mate `[s, e)` lies entirely within a single exon.
    fn holds(&self, s: i64, e: i64) -> bool {
        self.exons.iter().any(|&(es, ee)| es <= s && e <= ee)
    }
}

#[derive(Clone, Copy)]
struct Mate {
    start: i64,
    end: i64,
    tid: i32,
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

    let mut pending: HashMap<Vec<u8>, Mate> = HashMap::new();
    let mut rec = RawRecord::default();
    loop {
        if raw::read_record(reader.get_mut(), &mut rec)? == 0 {
            break;
        }
        let flags = rec.flags();
        if flags & SKIP != 0 || flags & PAIRED == 0 || rec.mapping_quality() < mapq {
            continue;
        }
        let tid = rec.reference_sequence_id();
        if tid < 0 {
            continue;
        }
        let start = i64::from(rec.alignment_start());
        let end = start + ref_span(rec.cigar_ops());
        let name = rec.name().to_vec();
        match pending.remove(&name) {
            None => {
                pending.insert(name, Mate { start, end, tid });
            }
            Some(mate) => {
                if mate.tid == tid {
                    record_fragment(
                        &transcripts,
                        &by_chrom,
                        &ref_names,
                        tid,
                        (start, end),
                        (mate.start, mate.end),
                        &mut sizes,
                    );
                }
            }
        }
    }
    Ok(FragmentSizes { transcripts, sizes })
}

fn record_fragment(
    transcripts: &[Transcript],
    by_chrom: &ChromIndex,
    ref_names: &[String],
    tid: i32,
    m1: (i64, i64),
    m2: (i64, i64),
    sizes: &mut [Vec<i64>],
) {
    let Some(chrom) = ref_names.get(tid as usize) else {
        return;
    };
    let Some((entries, max_span)) = by_chrom.get(chrom) else {
        return;
    };
    let left_start = m1.0.min(m2.0);
    let right_end = m1.1.max(m2.1);
    let hi = entries.partition_point(|&(s, _)| s <= left_start);
    for &(tx_start, idx) in entries[..hi].iter().rev() {
        if tx_start < left_start - *max_span {
            break;
        }
        let t = &transcripts[idx];
        if t.tx_end < right_end || !t.holds(m1.0, m1.1) || !t.holds(m2.0, m2.1) {
            continue;
        }
        // mRNA span plus one per exon junction the fragment crosses (RSeQC).
        if let (Some((a_off, a_ex)), Some((b_off, b_ex))) = (t.mrna(left_start), t.mrna(right_end))
        {
            sizes[idx].push((b_off - a_off) + (b_ex - a_ex) as i64);
        }
    }
}

fn ref_span<I: Iterator<Item = (u8, u32)>>(ops: I) -> i64 {
    let mut span = 0;
    for (op, len) in ops {
        if matches!(op, 0 | 2 | 3 | 7 | 8) {
            span += i64::from(len);
        }
    }
    span
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
