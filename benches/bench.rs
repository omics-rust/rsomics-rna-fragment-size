use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_frag(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-rna-fragment-size");
    let g = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    c.bench_function("rsomics-rna-fragment-size golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
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
        });
    });
}

criterion_group!(benches, bench_frag);
criterion_main!(benches);
