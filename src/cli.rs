use std::num::NonZero;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, Tool, ToolMeta};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-rna-fragment-size",
    version,
    about = "Per-transcript mRNA fragment-size distribution for paired RNA-seq — port of RSeQC RNA_fragment_size.py"
)]
pub struct Cli {
    /// Input BAM file.
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Reference gene model in 12-column BED format.
    #[arg(short = 'r', long = "refgene")]
    pub bed: PathBuf,

    /// Minimum mapping quality for an alignment to be "uniquely mapped".
    #[arg(long = "mapq", default_value_t = 30)]
    pub mapq: u8,

    /// Minimum number of fragments; transcripts below it report count then 0 0 0.
    #[arg(short = 'n', long = "frag-num", default_value_t = 3)]
    pub frag_num: usize,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let workers = self
            .common
            .threads
            .and_then(NonZero::new)
            .unwrap_or_else(|| {
                std::thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN)
            });
        let fs =
            rsomics_rna_fragment_size::fragment_sizes(&self.input, &self.bed, self.mapq, workers)?;
        let stdout = std::io::stdout();
        fs.write_report(stdout.lock(), self.frag_num)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
