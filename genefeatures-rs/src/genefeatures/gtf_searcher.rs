use crate::gtf_record::GtfRecord;
use std::collections::HashMap;

#[derive(Debug, PartialEq)]
pub enum GtfRecordValue <'a> {
    Str(&'a str),
    U32(u32),
    U8(u8),
    F64(f64),
    OptStr(Option<&'a str>),
    OptU8(Option<u8>),
    OptVecStr(Option<Vec<&'a str>>),
    NoMatch,
}


#[derive(Debug)]
pub struct GtfSearcher<'a> {
    condition: HashMap<&'a str, GtfRecordValue<'a>>,
}
impl <'a>GtfSearcher<'a> {
    pub fn new(condition: HashMap<&'a str, GtfRecordValue<'a>>) -> Self {
        GtfSearcher {
            condition,
        }
    }

    pub fn get_value_from_gtf_record(record: &'a GtfRecord, field: &str) -> GtfRecordValue<'a> {
        match field {
            "seqname" => GtfRecordValue::Str(&record.seqname),
            "source" => GtfRecordValue::Str(&record.source),
            "feature" => GtfRecordValue::Str(&record.feature),
            "start" => GtfRecordValue::U32(record.start),
            "end" => GtfRecordValue::U32(record.end),
            "score" => GtfRecordValue::F64(record.score),
            "frame" => GtfRecordValue::U8(record.frame),
            "gene_id" => GtfRecordValue::OptStr(record.gene_id.as_deref()),
            "gene_name" => GtfRecordValue::OptStr(record.gene_name.as_deref()),
            "gene_source" => GtfRecordValue::OptStr(record.gene_source.as_deref()),
            "gene_version" => GtfRecordValue::OptU8(record.gene_version),
            "gene_biotype" => GtfRecordValue::OptStr(record.gene_biotype.as_deref()),
            "transcript_id" => GtfRecordValue::OptStr(record.transcript_id.as_deref()),
            "transcript_name" => GtfRecordValue::OptStr(record.transcript_name.as_deref()),
            "transcript_source" => GtfRecordValue::OptStr(record.transcript_source.as_deref()),
            "transcript_version" => GtfRecordValue::OptU8(record.transcript_version),
            "transcript_biotype" => GtfRecordValue::OptStr(record.transcript_biotype.as_deref()),
            "transcript_support_level" => GtfRecordValue::OptU8(record.transcript_support_level),
            "exon_id" => GtfRecordValue::OptStr(record.exon_id.as_deref()),
            "exon_number" => GtfRecordValue::OptU8(record.exon_number),
            "exon_version" => GtfRecordValue::OptU8(record.exon_version),
            "protein_id" => GtfRecordValue::OptStr(record.protein_id.as_deref()),
            "protein_version" => GtfRecordValue::OptU8(record.protein_version),
            "ccds_id" => GtfRecordValue::OptStr(record.ccds_id.as_deref()),
            "tag" => GtfRecordValue::OptVecStr(
                record.tag
                    .as_ref()
                    .map(
                        |v|  v.iter().map(|s| &**s).collect()
                    )
                ),
            _ => GtfRecordValue::NoMatch,
        }
    }

    pub fn match_condition(
        condition: &GtfRecordValue,
        record_value: &GtfRecordValue,
    ) -> bool {
        match (condition, record_value) {
            (GtfRecordValue::Str(cond), GtfRecordValue::Str(val)) => cond == val,
            (GtfRecordValue::U32(cond), GtfRecordValue::U32(val)) => cond == val,
            (GtfRecordValue::U8(cond), GtfRecordValue::U8(val)) => cond == val,
            (GtfRecordValue::F64(cond), GtfRecordValue::F64(val)) => cond == val,
            (GtfRecordValue::OptStr(cond), GtfRecordValue::OptStr(val)) => cond == val,
            (GtfRecordValue::OptU8(cond), GtfRecordValue::OptU8(val)) => cond == val,
            (GtfRecordValue::OptVecStr(cond), GtfRecordValue::OptVecStr(val)) => {
                match (cond, val) {
                    (Some(c), Some(v)) => { 
                        c.iter().any(|e1| v.iter().any(|e2| e1 == e2))
                    },
                    _ => false
                } 
            },
            (_, GtfRecordValue::NoMatch) => false,
            _ => false,
        }
    }

    pub fn find_match(&mut self, record: &'a GtfRecord) -> bool {
        let mut all_match = true;
        for (field, condition) in &self.condition {
            let record_value = GtfSearcher::get_value_from_gtf_record(record, field);
            if !GtfSearcher::match_condition(condition, &record_value) {
                all_match = false;
                break;
            }
        }
        let result: bool = if all_match {
            true
        } else {
            false
        };
        result
    }
}
