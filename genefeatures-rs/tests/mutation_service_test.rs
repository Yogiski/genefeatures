extern crate genefeatures;
use bio_seq::prelude::*;
use genefeatures::mutation_service::{self, MutRegex};


#[cfg(test)]
mod tests {

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
        let mt: MutRegex = MutRegex::PntDel;
        assert!(mt.is_match("34_36del"));
        assert!(mt.is_match("34_36delGG"));
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
        assert_eq!(mt.get_recipe("34_36delGGTinsTAA"), ["34", "36", "GGT", "TAA"])
    }

    #[test]
    fn test_mutate_sequence_dup() {
    }
}