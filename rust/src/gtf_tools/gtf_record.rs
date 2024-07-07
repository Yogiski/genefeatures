#[derive(Debug)]
pub struct GtfRecord {
    pub seqname: Box<str>,
    pub source: Box<str>,
    pub feature: Box<str>,
    pub start: u32,
    pub end: u32,
    pub score: f64,
    pub strand: Box<str>,
    pub frame: u8,
    pub gene_id: Option<Box<str>>,
    pub gene_name: Option<Box<str>>,
    pub gene_source: Option<Box<str>>,
    pub gene_version: Option<u8>,
    pub gene_biotype: Option<Box<str>>,
    pub transcript_id: Option<Box<str>>,
    pub transcript_name: Option<Box<str>>,
    pub transcript_source: Option<Box<str>>,
    pub transcript_version: Option<u8>,
    pub transcript_biotype: Option<Box<str>>,
    pub transcript_support_level: Option<u8>,
    pub exon_id: Option<Box<str>>,
    pub exon_number: Option<u8>,
    pub exon_version: Option<u8>,
    pub protein_id: Option<Box<str>>,
    pub protein_version: Option<u8>,
    pub ccds_id: Option<Box<str>>,
    pub tag: Option<Vec<Box<str>>>
}

impl GtfRecord {

    pub fn new(
        seqname: Box<str>,
        source: Box<str>,
        feature: Box<str>,
        start: u32,
        end: u32,
        score: f64,
        strand: Box<str>,
        frame: u8,
        gene_id: Option<Box<str>>,
        gene_name: Option<Box<str>>,
        gene_source: Option<Box<str>>,
        gene_version: Option<u8>,
        gene_biotype: Option<Box<str>>,
        transcript_id: Option<Box<str>>,
        transcript_name: Option<Box<str>>,
        transcript_source: Option<Box<str>>,
        transcript_version: Option<u8>,
        transcript_biotype: Option<Box<str>>,
        transcript_support_level: Option<u8>,
        exon_id: Option<Box<str>>,
        exon_number: Option<u8>,
        exon_version: Option<u8>,
        protein_id: Option<Box<str>>,
        protein_version: Option<u8>,
        ccds_id: Option<Box<str>>,
        tag: Option<Vec<Box<str>>>
    ) -> Self {
        Self {
            seqname,
            source,
            feature,
            start,
            end,
            score,
            strand,
            frame,
            gene_id,
            gene_name,
            gene_source,
            gene_version,
            gene_biotype,
            transcript_id,
            transcript_name,
            transcript_source,
            transcript_version,
            transcript_biotype,
            transcript_support_level,
            exon_id,
            exon_number,
            exon_version,
            protein_id,
            protein_version,
            ccds_id,
            tag
        }
    }

    pub fn from_gtf_line(line: &str) -> Self {

        let fields: Vec<&str> = line.split("\t").collect();
        if fields.len() != 9 {
            panic!("Invalid GTF line: expected 9 fields, found {}", fields.len())
        }
        // get mandatory fields from line 
        let seqname: Box<str> = Box::from(fields[0]);
        let source: Box<str> = Box::from(fields[1]);
        let feature: Box<str> = Box::from(fields[2]);
        let start: u32 = fields[3]
            .parse::<u32>()
            .expect("Invalid start field");
        let end: u32 = fields[4]
            .parse::<u32>()
            .expect("Invalid end field");
        let score: f64 = match fields[5] {
            "." => 0.0,
            _ => fields[5]
                .parse::<f64>()
                .expect("Invalid score field")
        };
        let strand: Box<str> = Box::from(fields[6]);
        let frame: u8 = match fields[7] {
            "." => 0,
            _ => fields[7]
                .parse::<u8>()
                .expect("Invalid frame field")
        };
        // init optional fields as none
        let mut gene_id: Option<Box<str>> = None;
        let mut gene_name: Option<Box<str>> = None;
        let mut gene_source: Option<Box<str>> = None;
        let mut gene_version: Option<u8> = None;
        let mut gene_biotype: Option<Box<str>> = None;
        let mut transcript_id: Option<Box<str>> = None;
        let mut transcript_name: Option<Box<str>> = None;
        let mut transcript_source: Option<Box<str>> = None;
        let mut transcript_version: Option<u8> = None;
        let mut transcript_biotype: Option<Box<str>> = None;
        let mut transcript_support_level: Option<u8> = None;
        let mut exon_id: Option<Box<str>> = None;
        let mut exon_number: Option<u8> = None;
        let mut exon_version: Option<u8> = None;
        let mut protein_id: Option<Box<str>> = None;
        let mut protein_version: Option<u8> = None;
        let mut ccds_id: Option<Box<str>> = None;
        let mut tag: Option<Vec<Box<str>>> = None;

        // process attributes field to file optional fields
        let attributes: &str = fields[8];
        let attr_iter = attributes
            .split(";")
            .map(str::trim)
            .filter(|s| !s.is_empty());

        for attribute in  attr_iter {
            let mut iter: std::str::SplitN<char> = attribute.splitn(2, ' ');
            let key: &str = iter
                .next()
                .expect("Invalid attribute format");
            let value: Box<str> = iter
                .next()
                .expect("Invalid attribute format")
                .trim_matches('"')
                .into();

            match key {
                "gene_id" => gene_id = Some(value),
                "gene_name" => gene_name = Some(value),
                "gene_source" => gene_source = Some(value),
                "gene_version" => gene_version = value.parse::<u8>().ok(),
                "gene_biotype" => gene_biotype = Some(value),
                "transcript_id" => transcript_id = Some(value),
                "transcript_name" => transcript_name = Some(value),
                "transcript_source" => transcript_source = Some(value),
                "transcript_biotype" => transcript_biotype = Some(value),
                "transcript_version" => transcript_version = value.parse::<u8>().ok(),
                "transcript_support_level" => transcript_support_level = value.parse::<u8>().ok(),
                "exon_id" => exon_id = Some(value),
                "exon_number" => exon_number = value.parse::<u8>().ok(),
                "exon_version" => exon_version = value.parse::<u8>().ok(),
                "protein_id" => protein_id = Some(value),
                "protein_version" => protein_version = value.parse::<u8>().ok(),
                "ccds_id" => ccds_id = Some(value),
                "tag" => {
                    if let Some(mut tags) = tag {
                        tags.push(value);
                        tag = Some(tags);
                    } else {
                        tag = Some(vec![value]);
                    }
                }, 
                _ => {}
            }
        }
        GtfRecord::new(
            seqname, source, feature, start, end, score, strand, frame,
            gene_id, gene_name, gene_source, gene_version, gene_biotype,
            transcript_id, transcript_name, transcript_source, transcript_version,
            transcript_biotype, transcript_support_level,
            exon_id, exon_number, exon_version, protein_id, protein_version,
            ccds_id, tag
        )
    }
}