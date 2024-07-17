use clap::{ Parser, Subcommand};
use genefeatures::make_gf_fasta;

#[derive(Parser, Debug)]
#[command(name = "GeneFeatures")]
#[command(about = "A CLI tool for generating sets of sequences from GTF file annotations")]
struct Cli {
    // path to gtf file
    #[command(subcommand)]
    command: Commands,

}

#[derive(Subcommand, Debug)]
enum Commands {
    MakeGfFasta {
        #[arg(short, long, value_name = "GTF FILE", help = "Path to gtf file containing gene element annotations of fasta file")]
        gtf: String,
        // path to fasta file
        #[arg(short, long, value_name = "FASTA FILE", help = "Path to fasta file containg genome sequences")]
        fasta: String,
        // path to feature file
        #[arg(short = 'n', long, value_name = "Gene Transcript Mutations text file", help = "Path to 3 column tab separated file. col1 = Ensembl Gene Id, col2 = Ensembl Transcript ID, col3 = Mutation in dna change notation. Col1 required, Cols 2/3 optional.")]
        genes: String,
        // Model ID, name of output
        #[arg(short, long, value_name = "Model Id", help = "Name or ID of model. This will be used to name output files.")]
        model: Option<String>,
    },
    MakeCdsIndex {
        #[arg(short, long, value_name = "GTF FILE")]
        gtf: String,
    }
}

fn main() {
    let cli: Cli = Cli::parse();
    match &cli.command {
        Commands::MakeGfFasta { gtf, fasta, genes, model } => {
            make_gf_fasta::main(gtf, fasta, genes, model.as_ref().unwrap())
        }
        Commands::MakeCdsIndex { gtf } => {
            todo!()
        }
    }
}
