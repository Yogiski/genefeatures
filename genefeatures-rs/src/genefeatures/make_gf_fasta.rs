use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};
use crate::gtf_searcher::{GtfSearcher, GtfRecordValue};
use crate::gtf_tree::{GtfTree, Transcript, Node};


pub fn read_gene_file(genes: &Path) -> Vec<(Box<str>, Option<Box<str>>, Option<Box<str>>)> {

    fn check_empty_str(to_check: Option<&str>) -> Option<Box<str>> {
        match to_check {
            Some(s) if s.is_empty() => None,
            Some(s) => Some(s.to_string().into_boxed_str()),
            None => None
        }
    }

    let genes_file: File = File::open(genes)
        .expect("Failed to open genes file");
    let mut gene_tups: Vec<(Box<str>, Option<Box<str>>, Option<Box<str>>)> = Vec::new(); 

    for line in  BufReader::new(genes_file).lines() {
        let line: String = line.expect("failed to read genes line");
        let mut splt_iter: std::str::SplitN<&str> = line.splitn(3, ",");
        gene_tups.push(
            (
                splt_iter
                    .next()
                    .expect("gene column value missing!")
                    .to_string()
                    .into_boxed_str(),
                check_empty_str(splt_iter.next()),
                check_empty_str(splt_iter.next())
            )
        )
    }
    gene_tups
}

fn make_gtf_searchers<'a> (
    gene_tups: &'a Vec<(Box<str>, Option<Box<str>>, Option<Box<str>>)>
) -> Vec<GtfSearcher<'a>> {

    gene_tups
        .into_iter()
        .map(|(gid, tid, _)| {

            let mut condition: HashMap<&str, GtfRecordValue> = HashMap::new();
            condition.insert("feature", GtfRecordValue::Str("transcript"));
            condition.insert("gene_id", GtfRecordValue::OptStr(Some(gid)));
            match tid {
                Some(t) => condition.insert(
                    "transcript_id",
                    GtfRecordValue::OptStr(Some(t))
                ),
                None => condition.insert(
                    "tag",
                    GtfRecordValue::OptVecStr(
                        Some(vec![
                            "MANE_select",
                            "Ensembl_canonical"
                        ])
                    )
                )
            };
            GtfSearcher::new(condition)
        })
        .collect()
}


fn find_matching_records<'a>(gtf: &'a GtfTree, searchers: &mut Vec<GtfSearcher<'a>>) -> Vec<Option<&'a Transcript>> {
    searchers
        .iter_mut()
        .map(|s| gtf.find_transcript(s))
        .collect()
}


pub fn main(gtf: &String, fasta: &String, genes: &String, model: &String) {

    let fasta: &Path = Path::new(fasta);

    println!("Reading gene transcript mutations file: {}", genes);
    let genes:Vec<(Box<str>, Option<Box<str>>, Option<Box<str>>)>  = read_gene_file(Path::new(genes));
    let mut searchers: Vec<GtfSearcher> = make_gtf_searchers(&genes);

    println!("Reading gtf file: {}", gtf);
    let gtf:GtfTree  = GtfTree::parse_gtf_file(Path::new(gtf));
    let transcripts = find_matching_records(&gtf, &mut searchers);








}