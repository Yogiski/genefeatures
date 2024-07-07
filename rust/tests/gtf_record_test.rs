extern crate gtf_tools;
use gtf_tools::gtf_record::GtfRecord;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gtf_record_from_line() {

        let line: String = String::from(
            "1\thavana\tCDS\t3385149\t3385286\t.\t+\t0\tgene_id \"ENSG00000142611\"; gene_version \"17\"; transcript_id \"ENST00000514189\"; transcript_version \"5\"; exon_number \"4\"; gene_name \"PRDM16\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"PRDM16-208\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; protein_id \"ENSP00000421400\"; protein_version \"1\"; tag \"basic\"; tag \"GENCODE Primary\"; transcript_support_level \"5\";"
        ); 
        let record: GtfRecord = GtfRecord::from_gtf_line(&line);
        assert_eq!(record.seqname, Box::from("1"));
        assert_eq!(record.source, Box::from("havana"));
        assert_eq!(record.feature, Box::from("CDS"));
        assert_eq!(record.start, 3385149);
        assert_eq!(record.end, 3385286);
        assert_eq!(record.score, 0.0);
        assert_eq!(record.strand, Box::from("+"));
        assert_eq!(record.frame, 0);
        assert_eq!(record.gene_id.as_deref(), Some("ENSG00000142611"));
        assert_eq!(record.gene_version, Some(17));
        assert_eq!(record.transcript_id.as_deref(), Some("ENST00000514189"));
        assert_eq!(record.transcript_version, Some(5));
        assert_eq!(record.exon_number, Some(4));
        assert_eq!(record.gene_name.as_deref(), Some("PRDM16"));
        assert_eq!(record.gene_source.as_deref(), Some("ensembl_havana"));
        assert_eq!(record.gene_biotype.as_deref(), Some("protein_coding"));
        assert_eq!(record.transcript_name.as_deref(), Some("PRDM16-208"));
        assert_eq!(record.transcript_source.as_deref(), Some("havana"));
        assert_eq!(record.transcript_biotype.as_deref(), Some("protein_coding"));
        assert_eq!(record.protein_id.as_deref(), Some("ENSP00000421400"));
        assert_eq!(record.protein_version, Some(1));
        assert_eq!(record.tag, Some(vec![Box::from("basic"), Box::from("GENCODE Primary")]));
        assert_eq!(record.transcript_support_level, Some(5));
    }
}