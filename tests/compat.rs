//! Differential test against RSeQC `RNA_fragment_size.py`. `tests/golden/golden.txt`
//! is the exact stdout RSeQC 5.x produced for `reads.bam` + `genes.bed` (BED12 with
//! a 2-exon gene, a 3-exon gene exercising multi-junction fragments, and a gene with
//! no fragments; the BAM has spanning pairs plus a low-mapq pair that must be
//! dropped). We must reproduce it byte-for-byte (with `-n 1`).

use std::path::PathBuf;
use std::process::Command;

fn run(bam: &str, bed: &str, golden: &str) {
    let g = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let bin = env!("CARGO_BIN_EXE_rsomics-rna-fragment-size");
    let out = Command::new(bin)
        .args([
            "-i",
            g.join(bam).to_str().unwrap(),
            "-r",
            g.join(bed).to_str().unwrap(),
            "-n",
            "1",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let ours = String::from_utf8(out.stdout).unwrap();
    let want = std::fs::read_to_string(g.join(golden)).unwrap();
    assert_eq!(ours, want, "fragment-size report differs from RSeQC golden");
}

#[test]
fn matches_rseqc_rna_fragment_size() {
    run("reads.bam", "genes.bed", "golden.txt");
}

/// Crafted pairs exercise each way the size rule differs from a naive
/// fragment-span measurement, with expected values taken from RSeQC 5.0.4:
/// a mate whose alignment runs past `tx_end` (kept on `mate_start <= tx_end`),
/// an N-CIGAR read spanning an exon junction and a read pair buried in an
/// intron (both counted), soft-clips excluded and insertions included in
/// read1's query length, and a read1 lying downstream of its mate.
#[test]
fn matches_rseqc_crafted_edge_cases() {
    run("crafted.bam", "crafted.bed", "crafted.txt");
}
