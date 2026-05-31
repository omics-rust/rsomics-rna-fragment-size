# rsomics-rna-fragment-size

Per-transcript mRNA fragment-size distribution for paired-end RNA-seq — a fast
Rust port of RSeQC `RNA_fragment_size.py`.

```
rsomics-rna-fragment-size -i aligned.bam -r refgene.bed12 [--mapq 30] [-n 3]
```

For each transcript it pairs uniquely-mapped mates that both fall inside the
transcript's exons and measures the fragment length in transcript (mRNA)
coordinates — the genomic span with introns removed — then prints the count,
mean, median and population standard deviation. Transcripts with fewer than `-n`
fragments report `count 0 0 0`.

## Origin

Independent Rust reimplementation of RSeQC `RNA_fragment_size.py`, based on the
tool's documented behaviour, the SAM/BAM and BED12 format specs, and black-box
testing against the upstream binary (which fixed the mRNA fragment formula —
spliced span plus one per exon junction — and the pair filters). No upstream
source was read. Test fixtures are independently generated.

License: MIT OR Apache-2.0.
Upstream credit: RSeQC <https://rseqc.sourceforge.net/> (LGPL-2.1+).
