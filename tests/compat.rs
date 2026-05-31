//! Differential test against RSeQC `RNA_fragment_size.py`. `tests/golden/golden.txt`
//! is the exact stdout RSeQC 5.x produced for `reads.bam` + `genes.bed` (BED12 with
//! a 2-exon gene, a 3-exon gene exercising multi-junction fragments, and a gene with
//! no fragments; the BAM has spanning pairs plus a low-mapq pair that must be
//! dropped). We must reproduce it byte-for-byte (with `-n 1`).

use std::path::PathBuf;
use std::process::Command;

#[test]
fn matches_rseqc_rna_fragment_size() {
    let g = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let bin = env!("CARGO_BIN_EXE_rsomics-rna-fragment-size");
    let out = Command::new(bin)
        .args([
            "-i",
            g.join("reads.bam").to_str().unwrap(),
            "-r",
            g.join("genes.bed").to_str().unwrap(),
            "-n",
            "1",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let ours = String::from_utf8(out.stdout).unwrap();
    let golden = std::fs::read_to_string(g.join("golden.txt")).unwrap();
    assert_eq!(
        ours, golden,
        "fragment-size report differs from RSeQC golden"
    );
}
