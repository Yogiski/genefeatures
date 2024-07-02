#[derive(Debug)]
pub struct GtfRecord {
    pub seqname: String,
    pub source: String,
    pub feature: String,
    pub start: u32,
    pub end: u32,
    pub score: f64,
    pub strand: String,
    pub frame: u8,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub gene_source: Option<String>,
    pub gene_version: Option<u8>,
    pub gene_biotype: Option<String>,
    pub transcript_id: Option<String>,
    pub transcript_name: Option<String>,
    pub transcript_source: Option<String>,
    pub transcript_version: Option<u8>,
    pub transcript_biotype: Option<String>,
    pub transcript_support_level: Option<u8>,
    pub exon_id: Option<String>,
    pub exon_number: Option<u8>,
    pub exon_version: Option<u8>,
    pub protein_id: Option<String>,
    pub protein_version: Option<u8>,
    pub ccds_id: Option<String>,
    pub tag: Option<Vec<String>>
}

impl GtfRecord {
    pub fn new(
        seqname: String,
        source: String,
        feature: String,
        start: u32,
        end: u32,
        score: f64,
        strand: String,
        frame: u8,
        gene_id: Option<String>,
        gene_name: Option<String>,
        gene_source: Option<String>,
        gene_version: Option<u8>,
        gene_biotype: Option<String>,
        transcript_id: Option<String>,
        transcript_name: Option<String>,
        transcript_source: Option<String>,
        transcript_version: Option<u8>,
        transcript_biotype: Option<String>,
        transcript_support_level: Option<u8>,
        exon_id: Option<String>,
        exon_number: Option<u8>,
        exon_version: Option<u8>,
        protein_id: Option<String>,
        protein_version: Option<u8>,
        ccds_id: Option<String>,
        tag: Option<Vec<String>>
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
            panic!(
                "Invalid GTF line: expected 9 fields, found {}", fields.len()
            )
        }
        let seqname: String = fields[0].to_string();
        let source: String = fields[1].to_string();
        let feature: String = fields[2].to_string();
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
        let strand: String = fields[6].to_string();
        let frame: u8 = match fields[7] {
            "." => 0,
            _ => fields[7]
                .parse::<u8>()
                .expect("Invalid frame field")
        };
        let mut gene_id: Option<String> = None;
        let mut gene_name: Option<String> = None;
        let mut gene_source: Option<String> = None;
        let mut gene_version: Option<u8> = None;
        let mut gene_biotype: Option<String> = None;
        let mut transcript_id: Option<String> = None;
        let mut transcript_name: Option<String> = None;
        let mut transcript_source: Option<String> = None;
        let mut transcript_version: Option<u8> = None;
        let mut transcript_biotype: Option<String> = None;
        let mut transcript_support_level: Option<u8> = None;
        let mut exon_id: Option<String> = None;
        let mut exon_number: Option<u8> = None;
        let mut exon_version: Option<u8> = None;
        let mut protein_id: Option<String> = None;
        let mut protein_version: Option<u8> = None;
        let mut ccds_id: Option<String> = None;
        let mut tag: Option<Vec<String>> = None;

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
            let value: String = iter
                .next()
                .expect("Invalid attribute format")
                .trim_matches('"')
                .to_string();
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