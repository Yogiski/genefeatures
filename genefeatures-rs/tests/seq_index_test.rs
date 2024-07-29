extern crate genefeatures;
use genefeatures::gtf_record_value_enum::GtfRecordValue;
use genefeatures::gtf_searcher::GtfSearcher;
use genefeatures::gtf_tree::{GtfTree, Node, Transcript};
use genefeatures::seq_index::SeqIdx;
use bio_seq::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

#[cfg(test)]
mod tests {

    use super::*;

    fn setup_seq_index(strand: &str) -> (Transcript, SeqIdx) {

        let gtf: GtfTree = GtfTree::parse_gtf_file(Path::new("tests/data/hs_four_oncogenes.gtf"));
        let mut conditions:HashMap<&str, GtfRecordValue>  = HashMap::new();
        match strand {
            "+" => conditions.insert("gene_id", GtfRecordValue::OptStr(Some("ENSG00000146648"))),
            "-" => conditions.insert("gene_id", GtfRecordValue::OptStr(Some("ENSG00000133703"))),
            _ => panic!()
        };
        conditions.insert(
            "tag",
            GtfRecordValue::OptVecStr(Some(vec!["MANE_select", "Ensembl_canonical"]))
        );

        let mut searcher: GtfSearcher = GtfSearcher::new(conditions);
        let transcript: &Transcript = gtf
            .find_transcript(&mut searcher)
            .expect("no transcript_found");

        (transcript.clone(), SeqIdx::new(transcript))
    }

    fn read_forward_seq() -> Seq<Dna> {

        let seq_string: String = fs::read_to_string("tests/data/egfr_full_seq.txt")
            .unwrap()
            .lines()
            .map(str::trim)
            .collect::<String>();
        let seq: Seq<Dna> = seq_string.try_into().unwrap();
        seq
    }

    fn read_reverse_seq() -> Seq<Dna> {
        let seq_string: String = match fs::read_to_string("tests/data/kras_full_seq.txt") {
            Ok(s) => s,
            Err(_) => panic!("Failed to read tests/data/kras_full_seq.txt")
        };
        let seq: Seq<Dna> = seq_string.try_into().unwrap();
        seq.revcomp()
    }

    #[test]
    fn test_seq_index_forward_genomic_index() {
        let (transcript, seq_idx): (Transcript, SeqIdx) = setup_seq_index("+");
        let seq_len:u64 = transcript.record.end.abs_diff(transcript.record.start);
        assert_eq!(seq_idx.genomic_index.get(&transcript.record.start), Some(&0u64));
        assert_eq!(seq_idx.genomic_index.get(&transcript.record.end), Some(&seq_len));
    }

    #[test]
    fn test_seq_index_reverse_genomic_index() {
        let (transcript, seq_idx): (Transcript, SeqIdx) = setup_seq_index("-");
        let seq_len:u64 = transcript.record.end.abs_diff(transcript.record.start);
        assert_eq!(seq_idx.genomic_index.get(&transcript.record.end), Some(&0u64));
        assert_eq!(seq_idx.genomic_index.get(&transcript.record.start), Some(&seq_len));
    }

    #[test]
    fn test_forward_mut_geno_index_match() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("+");
        assert_eq!(seq_idx.genomic_index.values().min(), seq_idx.mutation_index.values().min());
        assert_eq!(seq_idx.genomic_index.values().max(), seq_idx.mutation_index.values().max());
    }

    #[test]
    fn test_reverse_mut_geno_index_match() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("-");
        assert_eq!(seq_idx.genomic_index.values().min(), seq_idx.mutation_index.values().min());
        assert_eq!(seq_idx.genomic_index.values().max(), seq_idx.mutation_index.values().max());
    }

    #[test]
    fn test_reverse_seq_idx_matches_seq_len() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("-");
        let seq: Seq<Dna> = read_reverse_seq();
        assert_eq!(seq_idx.mutation_index.len(), seq.len())
    }

    #[test]
    fn test_forward_seq_idx_matches_seq_len() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("+");
        let seq: Seq<Dna> = read_forward_seq();
        assert_eq!(seq_idx.mutation_index.len(), seq.len())
    }

    #[test]
    fn test_forward_seq_utr_tss_junction() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("+");
        assert_eq!(seq_idx.mutation_index["1"], seq_idx.mutation_index["-1"] + 1);
    }
    
    #[test]
    fn test_reverse_seq_utr_tss_junction() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("-");
        assert_eq!(seq_idx.mutation_index["1"], seq_idx.mutation_index["-1"] + 1);
    }

    #[test]
    fn test_forward_seq_tss_idx_seq() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("+");
        let seq: Seq<Dna> = read_forward_seq();
        let tss_i: usize = seq_idx.mutation_index["1"] as usize;
        let tss_f: usize = seq_idx.mutation_index["4"] as usize;
        let tss_seq: Seq<Dna> = dna!["ATG"].into();
        assert_eq!(seq[tss_i..tss_f], tss_seq)
    }

    #[test]
    fn test_reverse_tss_idx_seq() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("-");
        let seq: Seq<Dna>  = read_reverse_seq();
        let tss_i: usize = seq_idx.mutation_index["1"] as usize;
        let tss_f: usize = seq_idx.mutation_index["4"] as usize;
        let tss_seq: Seq<Dna> = dna!["ATG"].into();
        assert_eq!(seq[tss_i..tss_f], tss_seq)
    }

    #[test]
    fn test_reverse_stop_idx_seq() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("-");
        let seq: Seq<Dna>  = read_reverse_seq();
        let tss_i: usize = seq_idx.mutation_index["*1"]  as usize;
        let tss_f: usize = seq_idx.mutation_index["*4"] as usize;
        let tss_seq: Seq<Dna> = dna!["TAA"].into();
        assert_eq!(seq[tss_i..tss_f], tss_seq)
    }

    #[test]
    fn test_forward_stop_idx_seq() {
        let (_, seq_idx): (Transcript, SeqIdx) = setup_seq_index("+");
        let seq: Seq<Dna>  = read_forward_seq();
        let tss_i: usize = seq_idx.mutation_index["*1"]  as usize;
        let tss_f: usize = seq_idx.mutation_index["*4"] as usize;
        let tss_seq: Seq<Dna> = dna!["TGA"].into();
        assert_eq!(seq[tss_i..tss_f], tss_seq)
    }


}
