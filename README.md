# rsomics-rna-fragment-size

Per-transcript mRNA fragment-size distribution for paired-end RNA-seq — a fast
Rust port of RSeQC `RNA_fragment_size.py`.

```
rsomics-rna-fragment-size -i aligned.bam -r refgene.bed12 [--mapq 30] [-n 3]
```

For each transcript it measures one fragment per `read1` alignment (flag `0x40`)
that starts inside the transcript span and whose mate also starts inside it,
then prints the count, mean, median and population standard deviation.
Transcripts with fewer than `-n` fragments report `count 0 0 0`.

The fragment size reproduces RSeQC exactly:

```
size = |mrna(mate_start) - mrna(read_start)| + junctions_crossed + qlen(read1)
```

`mrna(g)` is the number of exonic bases before genomic position `g` (an intronic
position collapses to the start of the next exon); `junctions_crossed` is the
difference of the two positions' exon indices — RSeQC's output carries a
one-base-per-junction offset that this term matches; and `qlen(read1)` is read1's
query-aligned length, summing the `M/I/=/X` CIGAR ops (insertions counted,
soft-clips and deletions excluded). A read1 counts for a transcript when it is
uniquely mapped (mapq ≥ `--mapq`), its mate is mapped (a proper-pair flag is not
required), `tx_start ≤ read_start < tx_end`, and `tx_start ≤ mate_start ≤ tx_end`.
The mate's own alignment may run past the exons or past `tx_end`; only its start
position is used.

## Origin

Independent Rust reimplementation of RSeQC `RNA_fragment_size.py`, based on the
tool's documented behaviour, the SAM/BAM and BED12 format specs, and black-box
testing against the upstream binary, which fixed the size formula above and the
pair filters. No upstream source was read. Test fixtures are independently
generated.

License: MIT OR Apache-2.0.
Upstream credit: RSeQC <https://rseqc.sourceforge.net/> (LGPL-2.1+).
