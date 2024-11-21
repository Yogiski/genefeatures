extern crate genefeatures;
use bio_seq::prelude::*;
use genefeatures::gtf_record_value_enum::GtfRecordValue;
use genefeatures::gtf_searcher::GtfSearcher;
use genefeatures::gtf_tree::{GtfTree, Node, Transcript};
use genefeatures::seq_index::SeqIdx;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

pub fn setup_seq_index(strand: &str) -> (Transcript, SeqIdx) {
    let gtf: GtfTree = GtfTree::parse_gtf_file(Path::new("tests/data/hs_four_oncogenes.gtf"));
    let mut conditions: HashMap<&str, GtfRecordValue> = HashMap::new();
    match strand {
        "+" => conditions.insert("gene_id", GtfRecordValue::OptStr(Some("ENSG00000146648"))),
        "-" => conditions.insert("gene_id", GtfRecordValue::OptStr(Some("ENSG00000133703"))),
        _ => panic!(),
    };
    conditions.insert(
        "tag",
        GtfRecordValue::OptVecStr(Some(vec!["MANE_select", "Ensembl_canonical"])),
    );

    let mut searcher: GtfSearcher = GtfSearcher::new(conditions);
    let transcript: &Transcript = gtf
        .find_transcript(&mut searcher)
        .expect("no transcript_found");

    (transcript.clone(), SeqIdx::new(transcript))
}

pub fn read_forward_seq() -> Seq<Dna> {
    let seq_string: String = fs::read_to_string("tests/data/egfr_full_seq.txt")
        .unwrap()
        .lines()
        .map(str::trim)
        .collect::<String>();
    let seq: Seq<Dna> = seq_string.try_into().unwrap();
    seq
}

pub fn read_reverse_seq() -> Seq<Dna> {
    let seq_string: String = match fs::read_to_string("tests/data/kras_full_seq.txt") {
        Ok(s) => s,
        Err(_) => panic!("Failed to read tests/data/kras_full_seq.txt"),
    };
    let seq: Seq<Dna> = seq_string.try_into().unwrap();
    seq.revcomp()
}
