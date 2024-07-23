use crate::gtf_record::GtfRecord;


#[derive(Debug, PartialEq)]
pub enum GtfRecordValue <'a> {
    Str(&'a str),
    U64(u64),
    U8(u8),
    F64(f64),
    OptStr(Option<&'a str>),
    OptU8(Option<u8>),
    OptVecStr(Option<Vec<&'a str>>),
    NoMatch,
}

impl <'a>From<&'a str> for GtfRecordValue<'a> {
    fn from(val: &'a str) -> Self {
        GtfRecordValue::Str(val)
    }
} 
impl <'a>From<&'a Box<str>> for GtfRecordValue<'a> {
    fn from(val: &'a Box<str>) -> Self {
        GtfRecordValue::Str(val)
    }
} 
impl <'a>From<u64> for GtfRecordValue<'a> {
    fn from(val: u64) -> Self {
        GtfRecordValue::U64(val)
    }
} 
impl <'a>From<u8> for GtfRecordValue<'a> {
    fn from(val: u8) -> Self {
        GtfRecordValue::U8(val)
    }
} 
impl <'a>From<f64> for GtfRecordValue<'a> {
    fn from(val: f64) -> Self {
        GtfRecordValue::F64(val)
    }
}
impl <'a>From<Option<&'a str>> for GtfRecordValue<'a> {
    fn from(val: Option<&'a str>) -> Self {
        GtfRecordValue::OptStr(val)
    }
}
impl <'a>From<Option<u8>> for GtfRecordValue<'a> {
    fn from(val: Option<u8>) -> Self {
        GtfRecordValue::OptU8(val)
    }
}
impl <'a>From<Option<&'a Vec<Box<str>>>> for GtfRecordValue<'a> {
    fn from(val: Option<&'a Vec<Box<str>>>) -> Self {
        GtfRecordValue::OptVecStr( 
            val.map(|v| v.iter().map(|s| &**s).collect())
        )
    }
}

impl <'a>GtfRecordValue<'a> {

    pub fn get_value_from_gtf_record(record: &'a GtfRecord, field: &str) -> Self {

        match field {
            "seqname" => GtfRecordValue::from(&record.seqname),
            "source" => GtfRecordValue::from(&record.source),
            "feature" => GtfRecordValue::from(&record.feature),
            "start" => GtfRecordValue::from(record.start),
            "end" => GtfRecordValue::from(record.end),
            "score" => GtfRecordValue::from(record.score),
            "frame" => GtfRecordValue::from(record.frame),
            "gene_id" => GtfRecordValue::from(record.gene_id.as_deref()),
            "gene_name" => GtfRecordValue::from(record.gene_name.as_deref()),
            "gene_source" => GtfRecordValue::from(record.gene_source.as_deref()),
            "gene_version" => GtfRecordValue::from(record.gene_version),
            "gene_biotype" => GtfRecordValue::from(record.gene_biotype.as_deref()),
            "transcript_id" => GtfRecordValue::from(record.transcript_id.as_deref()),
            "transcript_name" => GtfRecordValue::from(record.transcript_name.as_deref()),
            "transcript_source" => GtfRecordValue::from(record.transcript_source.as_deref()),
            "transcript_version" => GtfRecordValue::from(record.transcript_version),
            "transcript_biotype" => GtfRecordValue::from(record.transcript_biotype.as_deref()),
            "transcript_support_level" => GtfRecordValue::from(record.transcript_support_level),
            "exon_id" => GtfRecordValue::from(record.exon_id.as_deref()),
            "exon_number" => GtfRecordValue::from(record.exon_number),
            "exon_version" => GtfRecordValue::from(record.exon_version),
            "protein_id" => GtfRecordValue::from(record.protein_id.as_deref()),
            "protein_version" => GtfRecordValue::from(record.protein_version),
            "ccds_id" => GtfRecordValue::from(record.ccds_id.as_deref()),
            "tag" => GtfRecordValue::from(record.tag.as_ref()),
            _ => GtfRecordValue::NoMatch,
        }
    }
}