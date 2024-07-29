use crate::gtf_record::GtfRecord;
use crate::gtf_record_value_enum::GtfRecordValue;
use std::collections::HashMap;

#[derive(Debug)]
pub struct GtfSearcher<'a> {
    condition: HashMap<&'a str, GtfRecordValue<'a>>,
}
impl<'a> GtfSearcher<'a> {
    pub fn new(condition: HashMap<&'a str, GtfRecordValue<'a>>) -> Self {
        GtfSearcher { condition }
    }

    pub fn match_condition(condition: &GtfRecordValue, record_value: &GtfRecordValue) -> bool {
        match (condition, record_value) {
            (GtfRecordValue::Str(cond), GtfRecordValue::Str(val)) => cond == val,
            (GtfRecordValue::U64(cond), GtfRecordValue::U64(val)) => cond == val,
            (GtfRecordValue::U8(cond), GtfRecordValue::U8(val)) => cond == val,
            (GtfRecordValue::F64(cond), GtfRecordValue::F64(val)) => cond == val,
            (GtfRecordValue::OptStr(cond), GtfRecordValue::OptStr(val)) => cond == val,
            (GtfRecordValue::OptU8(cond), GtfRecordValue::OptU8(val)) => cond == val,
            (GtfRecordValue::OptVecStr(cond), GtfRecordValue::OptVecStr(val)) => {
                match (cond, val) {
                    (Some(c), Some(v)) => c.iter().any(|e1| v.iter().any(|e2| e1 == e2)),
                    _ => false,
                }
            }
            (_, GtfRecordValue::NoMatch) => false,
            _ => false,
        }
    }

    pub fn find_match(&self, record: &'a GtfRecord) -> bool {
        let mut all_match: bool = true;
        for (field, condition) in &self.condition {
            let record_value: GtfRecordValue =
                GtfRecordValue::get_value_from_gtf_record(record, field);
            if !GtfSearcher::match_condition(condition, &record_value) {
                all_match = false;
                break;
            }
        }
        all_match
    }
}
