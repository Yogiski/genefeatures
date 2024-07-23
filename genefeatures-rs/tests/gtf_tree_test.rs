use std::collections::HashMap;
use std::path::Path;

extern crate genefeatures;
use genefeatures::gtf_record_value_enum::GtfRecordValue;
use genefeatures::gtf_searcher::GtfSearcher;
use genefeatures::gtf_tree::{Cds, Contig, Gene, GtfTree, Node, NonCds, Transcript};

#[cfg(test)]
mod tests {

    use super::*;

    fn read_gtf_file() -> GtfTree {
        GtfTree::parse_gtf_file(Path::new("tests/data/hs_four_oncogenes.gtf"))
    }

    #[test]
    fn test_gtf_tree_from_file() {
        let gtf: GtfTree = read_gtf_file();
        assert!(!gtf.genebuild_last_updated.is_none());
        assert!(!gtf.genome_build.is_none());
        assert!(!gtf.contigs.is_empty());
        assert!(!gtf.contigs[0].genes.is_empty());
        assert!(!gtf.contigs[0].genes[0].transcripts.is_empty());
        assert!(!gtf.contigs[0].genes[0].transcripts[0]
            .cds_records
            .records
            .is_empty());
    }

    #[test]
    fn test_gtf_tree_find_transcript() {
        let gtf: GtfTree = read_gtf_file();
        let mapping: HashMap<&str, GtfRecordValue> = vec![
            ("gene_id", GtfRecordValue::OptStr(Some("ENSG00000133703"))),
            (
                "transcript_id",
                GtfRecordValue::OptStr(Some("ENST00000311936")),
            ),
        ]
        .into_iter()
        .collect();
        let mut searcher: GtfSearcher = GtfSearcher::new(mapping);
        let search_res: Option<&Transcript> = gtf.find_transcript(&mut searcher);
        assert!(search_res.is_some());

        let transcript: &Transcript = search_res.unwrap();
        assert!(transcript.record.gene_name.is_some());
        assert_eq!(
            search_res.unwrap().record.gene_name,
            Some("KRAS".to_string().into_boxed_str())
        );
        assert!(!transcript.cds_records.records.is_empty())
    }
}
