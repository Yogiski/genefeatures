# genefeatures

Primary functionality of GeneFeatures is to be able to generate fasta files with mutated dna sequences. GeneFeatures takes in a set of genes/transcripts and mutations with a corresponding reference fasta and gtf and returns a fasta of mutated dna sequences. Mutation formats accepted are DNA change and protein change, respective examples: g.76C>T, p.G12C.

Currently there is no tool to generate mutated consensus or reference sequences when the user does not have access to bam files or vcf files generated from a whole genome or whole exome sequencing experiment. GeneFeatures fills this niche allowing for users to apply and analyze mutated reference sequence to match the state of a given model of interest. 

Genefeatures is a work in progress and does not achieve the complete functionality described above.
It is currently being refactored from python to rust for performance improvements.

If this seems like an interesting or useful tool kindly reach out to my email: cyogodzi@gmail, with any questions or comments you may have. Or submit an issue with any features you'd like to see implemented.

 