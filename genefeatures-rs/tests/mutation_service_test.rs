extern crate genefeatures;
use genefeatures::gtf_tree::Transcript;
use genefeatures::mutation_service::{self, MutRegex};
use genefeatures::seq_index::SeqIdx;
mod test_utils;
use bio_seq::prelude::*;

#[cfg(test)]
mod tests {

    use mutation_service::process_mutation;

    use super::*;

    #[test]
    fn test_sub_regex_match() {
        let mt: MutRegex = MutRegex::Sub;
        assert!(mt.is_match("76A>T"));
    }

    #[test]
    fn test_sub_regex_get_recipe() {
        let mt: MutRegex = MutRegex::Sub;
        assert_eq!(mt.get_recipe("76A>T"), vec!["76", "A", "T"]);
    }

    #[test]
    fn test_pnt_del_regex_match() {
        let mt: MutRegex = MutRegex::PntDel;
        assert!(mt.is_match("34del"));
        assert!(mt.is_match("34delG"));
        assert!(!mt.is_match("31_34delG"));
    }

    #[test]
    fn test_pnt_del_regext_get_recipe() {
        let mt: MutRegex = MutRegex::PntDel;
        let recipe: Vec<&str> = mt.get_recipe("34del");
        assert_eq!(recipe, ["34", ""]);
        let recipe: Vec<&str> = mt.get_recipe("34delG");
        assert_eq!(recipe, ["34", "G"]);
    }

    #[test]
    fn test_rng_del_match() {
        let mt: MutRegex = MutRegex::RngDel;
        assert!(mt.is_match("34_36del"));
        assert!(mt.is_match("34_36delGG"));
        assert!(!mt.is_match("36del"));
    }

    #[test]
    fn test_rng_del_regex_get_recipe() {
        let mt: MutRegex = MutRegex::RngDel;
        let recipe: Vec<&str> = mt.get_recipe("34_36del");
        assert_eq!(recipe, ["34", "36", ""]);
        let recipe: Vec<&str> = mt.get_recipe("34_36delGG");
        assert_eq!(recipe, ["34", "36", "GG"]);
    }

    #[test]
    fn test_ins_regex_match() {
        let mt: MutRegex = MutRegex::Ins;
        assert!(mt.is_match("34_35insAAA"));
    }
    #[test]
    fn test_ins_regex_get_recipe() {
        let mt: MutRegex = MutRegex::Ins;
        assert_eq!(mt.get_recipe("34_35insAAA"), ["34", "35", "AAA"]);
    }

    #[test]
    fn test_dup_regex_match() {
        let mt: MutRegex = MutRegex::Dup;
        assert!(mt.is_match("34_36dup"));
        assert!(mt.is_match("34_36dupGGT"));
    }

    #[test]
    fn test_dup_regex_get_recipe() {
        let mt: MutRegex = MutRegex::Dup;
        assert_eq!(mt.get_recipe("34_36dup"), ["34", "36", ""]);
        assert_eq!(mt.get_recipe("34_36dupGGT"), ["34", "36", "GGT"])
    }

    #[test]
    fn test_inv_regex_match() {
        let mt: MutRegex = MutRegex::Inv;
        assert!(mt.is_match("34_36inv"));
        assert!(mt.is_match("34_36inv3"));
    }
    #[test]
    fn test_inv_regex_get_recipe() {
        let mt: MutRegex = MutRegex::Inv;
        assert_eq!(mt.get_recipe("34_36inv"), ["34", "36", ""]);
        assert_eq!(mt.get_recipe("34_36inv3"), ["34", "36", "3"]);
    }

    #[test]
    fn test_indel_regex_match() {
        let mt: MutRegex = MutRegex::InDel;
        assert!(mt.is_match("34_36delinsTAA"));
        assert!(mt.is_match("34_36delGGTinsTAA"))
    }
    #[test]
    fn test_index_regex_match() {
        let mt: MutRegex = MutRegex::InDel;
        assert_eq!(mt.get_recipe("34_36delinsTAA"), ["34", "36", "TAA"]);
        assert_eq!(
            mt.get_recipe("34_36delGGTinsTAA"),
            ["34", "36", "GGT", "TAA"]
        )
    }

    #[test]
    fn test_mutate_sequence_sub_forward() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("+");
        let seq: Seq<Dna> = test_utils::read_forward_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("3G>T".to_string(), seq, seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("ATT").into();
        assert_ne!(wt_start, mt_start);
        assert_eq!(target, mt_start)
    }

    #[test]
    fn test_mutate_sequence_sub_reverse() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("-");
        let seq: Seq<Dna> = test_utils::read_reverse_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("3G>T".to_string(), seq, seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("ATT").into();
        assert_ne!(wt_start, mt_start);
        assert_eq!(target, mt_start)
    }

    #[test]
    fn test_mutate_sequence_pnt_del_forward() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("+");
        let seq: Seq<Dna> = test_utils::read_forward_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("2del".to_string(), seq.clone(), seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("AGC").into();

        assert_ne!(seq.len(), mt_seq.len().clone());
        assert_eq!(mt_seq.len(), seq.len() - 1);
        assert_ne!(wt_start.to_string(), mt_start.to_string());
        assert_eq!(mt_start.to_string(), target.to_string());
    }

    #[test]
    fn test_mutate_sequence_pnt_del_reverse() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("-");
        let seq: Seq<Dna> = test_utils::read_reverse_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("2del".to_string(), seq.clone(), seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("AGA").into();

        assert_ne!(seq.len(), mt_seq.len().clone());
        assert_eq!(mt_seq.len(), seq.len() - 1);
        assert_ne!(wt_start.to_string(), mt_start.to_string());
        assert_eq!(mt_start.to_string(), target.to_string());
    }

    #[test]
    fn test_mutate_sequence_range_del_forward() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("+");
        let seq: Seq<Dna> = test_utils::read_forward_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("1_3del".to_string(), seq.clone(), seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("CGA").into();

        assert_ne!(seq.len(), mt_seq.len().clone());
        assert_eq!(mt_seq.len(), seq.len() - 3);
        assert_ne!(wt_start.to_string(), mt_start.to_string());
        assert_eq!(mt_start.to_string(), target.to_string());
    }

    #[test]
    fn test_mutate_sequence_range_del_reverse() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("-");
        let seq: Seq<Dna> = test_utils::read_reverse_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("1_3del".to_string(), seq.clone(), seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("ACT").into();

        assert_ne!(seq.len(), mt_seq.len().clone());
        assert_eq!(mt_seq.len(), seq.len() - 3);
        assert_ne!(wt_start.to_string(), mt_start.to_string());
        assert_eq!(mt_start.to_string(), target.to_string());
    }

    #[test]
    fn test_mutate_sequence_dup_forward() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("+");
        let seq: Seq<Dna> = test_utils::read_forward_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let dup_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["4"] as usize)..=(seq_idx.mutation_index["6"] as usize);

        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let wt_dup: Seq<Dna> = seq[dup_codon.clone()].to_owned();

        let mt_seq: Seq<Dna> = process_mutation("1_3dup".to_string(), seq, seq_idx);

        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let mt_dup: Seq<Dna> = mt_seq[dup_codon.clone()].to_owned();

        assert_eq!(wt_start, mt_start);
        assert_ne!(wt_dup, mt_dup);
        assert_eq!(wt_start, mt_dup);
    }

    #[test]
    fn test_mutate_sequence_dup_reverse() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("-");
        let seq: Seq<Dna> = test_utils::read_reverse_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let dup_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["4"] as usize)..=(seq_idx.mutation_index["6"] as usize);

        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let wt_dup: Seq<Dna> = seq[dup_codon.clone()].to_owned();

        let mt_seq: Seq<Dna> = process_mutation("1_3dup".to_string(), seq, seq_idx);

        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let mt_dup: Seq<Dna> = mt_seq[dup_codon.clone()].to_owned();

        assert_eq!(wt_start, mt_start);
        assert_ne!(wt_dup, mt_dup);
        assert_eq!(wt_start, mt_dup);
    }

    #[test]
    fn test_mutate_sequence_ins_forward() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("+");
        let seq: Seq<Dna> = test_utils::read_forward_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let dup_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["4"] as usize)..=(seq_idx.mutation_index["6"] as usize);

        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let wt_ins: Seq<Dna> = seq[dup_codon.clone()].to_owned();

        let mt_seq: Seq<Dna> = process_mutation("3_4insATG".to_string(), seq, seq_idx);

        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let mt_ins: Seq<Dna> = mt_seq[dup_codon.clone()].to_owned();

        assert_eq!(wt_start, mt_start);
        assert_ne!(wt_ins, mt_ins);
        assert_eq!(wt_start, mt_ins);
    }

    #[test]
    fn test_mutate_sequence_ins_reverse() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("-");
        let seq: Seq<Dna> = test_utils::read_reverse_seq();

        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let dup_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["4"] as usize)..=(seq_idx.mutation_index["6"] as usize);

        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let wt_ins: Seq<Dna> = seq[dup_codon.clone()].to_owned();

        let mt_seq: Seq<Dna> = process_mutation("3_4insATG".to_string(), seq, seq_idx);

        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let mt_ins: Seq<Dna> = mt_seq[dup_codon.clone()].to_owned();

        assert_eq!(wt_start, mt_start);
        assert_ne!(wt_ins, mt_ins);
        assert_eq!(wt_start, mt_ins);
    }

    #[test]
    fn test_mutate_sequence_inv_forward() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("+");
        let seq: Seq<Dna> = test_utils::read_forward_seq();
        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("1_3inv3".to_string(), seq.clone(), seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("GTA").into();

        assert_ne!(wt_start, mt_start);
        assert_eq!(mt_start, target);
    }

    #[test]
    fn test_mutate_sequence_inv_reverse() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("-");
        let seq: Seq<Dna> = test_utils::read_reverse_seq();
        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["3"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("1_3inv3".to_string(), seq.clone(), seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("GTA").into();

        assert_ne!(wt_start, mt_start);
        assert_eq!(mt_start, target);
    }

    #[test]
    fn test_mutate_sequence_indel_forward() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("+");
        let seq: Seq<Dna> = test_utils::read_forward_seq();
        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["6"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("1_3delinsAAACCC".to_string(), seq.clone(), seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("AAACCC").into();

        assert_eq!(seq.len() + 3, mt_seq.len());
        assert_ne!(wt_start, mt_start);
        assert_eq!(mt_start, target);
    }

    #[test]
    fn test_mutate_sequence_indel_reverse() {
        let (_, seq_idx): (Transcript, SeqIdx) = test_utils::setup_seq_index("-");
        let seq: Seq<Dna> = test_utils::read_reverse_seq();
        let start_codon: std::ops::RangeInclusive<usize> =
            (seq_idx.mutation_index["1"] as usize)..=(seq_idx.mutation_index["6"] as usize);
        let wt_start: Seq<Dna> = seq[start_codon.clone()].to_owned();
        let mt_seq: Seq<Dna> = process_mutation("1_3delinsAAACCC".to_string(), seq.clone(), seq_idx);
        let mt_start: Seq<Dna> = mt_seq[start_codon.clone()].to_owned();
        let target: Seq<Dna> = dna!("AAACCC").into();

        assert_eq!(seq.len() + 3, mt_seq.len());
        assert_ne!(wt_start, mt_start);
        assert_eq!(mt_start, target);
    }
}
