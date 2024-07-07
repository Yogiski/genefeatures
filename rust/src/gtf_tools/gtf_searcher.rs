use std::collections::HashMap;
use gtf_record::GtfRecord;

#[derive(Debug, PartialEq)]
pub enum GtfRecordValue<'a> {
    BoxStr(&'a Box<str>),
    String(String),
    U32(u32),
    U8(u8),
    F64(f64),
    VecString(Vec<Box<str>>),
    OptBoxStr(&'a Option<Box<str>>),
    OptU8(Option<u8>),
    OptVecBoxStr(&'a Option<Vec<Box<str>>>),
    NoMatch
}
impl <'a>GtfRecordValue<'a> {

    fn eq_option_box_str(&self, compare: &Option<Box<str>>) -> bool {
        match self {
            GtfRecordValue::OptBoxStr(Some(val)) => compare.as_deref() == Some(val.as_ref()),
            GtfRecordValue::OptBoxStr(None) => compare.is_none(),
            _ => false
        }
    }
    fn eq_option_u8(& self, compare: &Option<u8>) -> bool {
        match self {
            GtfRecordValue::OptU8(val) => compare == val,
            _ => false
        }
    }
    fn eq_option_vec_box(&self, compare: &Option<Vec<Box<str>>>) -> bool {
        match self {
            GtfRecordValue::OptVecBoxStr(Some(val)) => compare.as_ref() == Some(val),
            GtfRecordValue::OptVecBoxStr(None) => compare.is_none(),
            _ => false
        }
    }
}

#[derive(Debug)]
pub struct GtfSearcher<'a> {
    condition: HashMap<&'a str, GtfRecordValue<'a>>,
    matches: Vec<&'a GtfRecord>
}
impl <'a>GtfSearcher<'a> {

    fn new(condition: HashMap<&'a str, GtfRecordValue<'a>>) -> Self {
        GtfSearcher {
            condition,
            matches: Vec::new()
        }
    }

    fn get_value_from_gtf_record(record: &'a GtfRecord, field: &str) -> GtfRecordValue<'a> {
        match field {
            "seqname" => GtfRecordValue::BoxStr(&record.seqname),
            "source" => GtfRecordValue::BoxStr(&record.seqname),
            "feature" => GtfRecordValue::BoxStr(&record.feature),
            "start" => GtfRecordValue::U32(record.start),
            "end" => GtfRecordValue::U32(record.end),
            "score" => GtfRecordValue::F64(record.score),
            "frame" => GtfRecordValue::U8(record.frame),
            "gene_id" => GtfRecordValue::OptBoxStr(&record.gene_id),
            "gene_name" => GtfRecordValue::OptBoxStr(&record.gene_name),
            "gene_source" => GtfRecordValue::OptBoxStr(&record.gene_source),
            "gene_version" => GtfRecordValue::OptU8(record.gene_version),
            "gene_biotype" => GtfRecordValue::OptBoxStr(&record.gene_biotype),
            "transcript_id" => GtfRecordValue::OptBoxStr(&record.transcript_id),
            "transcript_name" => GtfRecordValue::OptBoxStr(&record.transcript_name),
            "transcript_source" => GtfRecordValue::OptBoxStr(&record.transcript_source),
            "transcript_version" => GtfRecordValue::OptU8(record.transcript_version),
            "transcript_biotype" => GtfRecordValue::OptBoxStr(&record.transcript_biotype),
            "transcript_support_level" => GtfRecordValue::OptU8(record.transcript_support_level),
            "exon_id" => GtfRecordValue::OptBoxStr(&record.exon_id),
            "exon_number" => GtfRecordValue::OptU8(record.exon_number),
            "exon_version" => GtfRecordValue::OptU8(record.exon_version),
            "protein_id" => GtfRecordValue::OptBoxStr(&record.protein_id),
            "protein_version" => GtfRecordValue::OptU8(record.protein_version),
            "ccds_id" => GtfRecordValue::OptBoxStr(&record.ccds_id),
            "tag" => GtfRecordValue::OptVecBoxStr(&record.tag),
            _ => GtfRecordValue::NoMatch
        }

    }

    fn match_condition(condition: &GtfRecordValue<'a>, record_value: &GtfRecordValue<'a>) -> bool {
        match (condition, record_value) {
            (GtfRecordValue::BoxStr(cond), GtfRecordValue::BoxStr(val)) => cond == val,
            (GtfRecordValue::U32(cond), GtfRecordValue::U32(val)) => cond == val,
            (GtfRecordValue::U8(cond), GtfRecordValue::U8(val)) => cond == val,
            (GtfRecordValue::F64(cond), GtfRecordValue::F64(val)) => cond == val,
            (GtfRecordValue::OptBoxStr(cond), GtfRecordValue::OptBoxStr(val)) =>cond == val,
            (GtfRecordValue::OptU8(cond), GtfRecordValue::OptU8(val)) => cond == val,
            (GtfRecordValue::OptVecBoxStr(cond), GtfRecordValue::OptVecBoxStr(val)) => cond == val,
            (_, GtfRecordValue::NoMatch) => false,
            _ => false
        }
    }

    fn find_match(&mut self, record: &'a GtfRecord) -> Vec<&'a GtfRecord> {
        let mut all_match = true;
        for (field, condition) in &self.condition {
            let record_value = GtfSearcher::get_value_from_gtf_record(record, field);
            if !GtfSearcher::match_condition(condition, &record_value) {
                all_match = false;
                break;
            }
        }
        if all_match {
            self.matches.push(record);
        }
        self.matches.clone()
    }
}