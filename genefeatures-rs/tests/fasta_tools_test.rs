extern crate genefeatures;
use genefeatures::fasta_tools::call_samtools_faidx;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_call_samtools_faidx() {

        let res_seq = call_samtools_faidx(
            "tests/data/trunc_hs.grch38.dna.chr1.fa", "1", 3069260, 3069262
        );

        match res_seq {
            Ok(seq) => assert_eq!(seq, "ATG"),
            Err(e) => {
                eprintln!("error extracting sequence {}", e);
                assert!(false)
            } 
        }
    }
}