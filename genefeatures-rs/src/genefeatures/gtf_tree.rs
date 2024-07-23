use crate::gtf_record::GtfRecord;
use crate::gtf_searcher::GtfSearcher;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

pub trait Node {
    fn add_record(&mut self, _record: GtfRecord) {}
    fn process_staged_records(&mut self) {}
    fn find_transcript<'a>(&'a self, _searcher: &mut GtfSearcher<'a>) -> Option<&Transcript> {
        None
    }
}

#[derive(Debug)]
pub struct Cds {
    pub records: Vec<GtfRecord>,
}
impl Cds {
    fn new() -> Self {
        Cds {
            records: Vec::new(),
        }
    }
}
impl Node for Cds {
    fn add_record(&mut self, record: GtfRecord) {
        self.records.push(record)
    }
}

#[derive(Debug)]
pub struct NonCds {
    pub records: Vec<GtfRecord>,
}
impl NonCds {
    fn new() -> Self {
        NonCds {
            records: Vec::new(),
        }
    }
}
impl Node for NonCds {
    fn add_record(&mut self, record: GtfRecord) {
        self.records.push(record)
    }
}

#[derive(Debug)]
pub struct Transcript {
    pub record: GtfRecord,
    pub cds_records: Cds,
    pub non_cds_records: NonCds,
}
impl Transcript {
    fn new(record: GtfRecord) -> Self {
        Transcript {
            record,
            cds_records: Cds::new(),
            non_cds_records: NonCds::new(),
        }
    }
}
impl Node for Transcript {
    fn add_record(&mut self, record: GtfRecord) {
        match record.feature.as_ref() {
            "CDS" => self.cds_records.add_record(record),
            _ => self.non_cds_records.add_record(record),
        }
    }
}

#[derive(Debug)]
pub struct Gene {
    pub record: GtfRecord,
    pub transcripts: Vec<Transcript>,
    staging: VecDeque<GtfRecord>,
}
impl Gene {
    fn new(record: GtfRecord) -> Self {
        Gene {
            record,
            transcripts: Vec::new(),
            staging: VecDeque::new(),
        }
    }
}
impl Node for Gene {
    fn add_record(&mut self, record: GtfRecord) {
        if let Some(transcript) = self
            .transcripts
            .iter_mut()
            .find(|t| t.record.transcript_id == record.transcript_id)
        {
            transcript.add_record(record)
        } else {
            match record.feature.as_ref() {
                "transcript" => self.transcripts.push(Transcript::new(record)),
                _ => self.staging.push_back(record),
            }
        }
    }
    fn process_staged_records(&mut self) {
        while let Some(record) = self.staging.pop_front() {
            self.add_record(record);
        }
    }

    fn find_transcript<'a>(&'a self, searcher: &mut GtfSearcher<'a>) -> Option<&Transcript> {
        self.transcripts
            .iter()
            .find(|t| searcher.find_match(&t.record))
    }
}

#[derive(Debug)]
pub struct Contig {
    pub genes: Vec<Gene>,
    pub name: String,
    staging: VecDeque<GtfRecord>,
}
impl Contig {
    fn new(name: String) -> Self {
        Contig {
            name,
            genes: Vec::new(),
            staging: VecDeque::new(),
        }
    }
}
impl Node for Contig {
    fn add_record(&mut self, record: GtfRecord) {
        if let Some(gene) = self
            .genes
            .iter_mut()
            .find(|g| g.record.gene_id == record.gene_id)
        {
            gene.add_record(record)
        } else {
            match record.feature.as_ref() {
                "gene" => self.genes.push(Gene::new(record)),
                _ => self.staging.push_back(record),
            }
        }
    }
    fn process_staged_records(&mut self) {
        while let Some(record) = self.staging.pop_front() {
            self.add_record(record);
        }
    }
    fn find_transcript<'a>(&'a self, searcher: &mut GtfSearcher<'a>) -> Option<&Transcript> {
        self.genes.iter().find_map(|g| g.find_transcript(searcher))
    }
}

#[derive(Debug)]
pub struct GtfTree {
    pub contigs: Vec<Contig>,
    pub genome_build: Option<String>,
    pub genome_version: Option<String>,
    pub genome_date: Option<String>,
    pub genome_build_accension: Option<String>,
    pub genebuild_last_updated: Option<String>,
}
impl GtfTree {
    pub fn new() -> Self {
        GtfTree {
            contigs: Vec::new(),
            genome_build: None,
            genome_version: None,
            genome_date: None,
            genome_build_accension: None,
            genebuild_last_updated: None,
        }
    }
    pub fn add_contig(&mut self, contig: Contig) {
        self.contigs.push(contig)
    }
    pub fn add_metadata(&mut self, key: &str, val: &str) {
        match key {
            "genome-build" => self.genome_build = Some(val.to_string()),
            "genome-version" => self.genome_version = Some(val.to_string()),
            "genome-date" => self.genome_date = Some(val.to_string()),
            "genome-build-accession" => self.genome_build_accension = Some(val.to_string()),
            "genebuild-last-updated" => self.genebuild_last_updated = Some(val.to_string()),
            _ => {}
        }
    }
    pub fn parse_gtf_file(gtf_path: &Path) -> Self {
        let gtf: File = File::open(gtf_path).expect("Failed to open GTF file");
        let mut gtf_tree: GtfTree = GtfTree::new();

        for line in io::BufReader::new(gtf).lines() {
            let to_record = line.expect("Failed to read line");
            if to_record.starts_with("#!") {
                let mut iter: std::str::SplitN<&str> = to_record
                    .strip_prefix("#!")
                    .expect("invalid metadata line!")
                    .splitn(2, " ");
                let key: &str = iter.next().expect("invalid metadata key");
                let val: &str = iter.next().expect("invalid metadata value");
                gtf_tree.add_metadata(key, val);
            } else {
                let record = GtfRecord::from_gtf_line(&to_record);
                gtf_tree.add_record(record)
            }
        }
        gtf_tree.process_staged_records();
        gtf_tree
    }
}
impl Node for GtfTree {
    fn add_record(&mut self, record: GtfRecord) {
        let seqname: String = record.seqname.clone().into_string();
        if let Some(contig) = self.contigs.iter_mut().find(|c| c.name == seqname) {
            contig.add_record(record)
        } else {
            let mut new_contig = Contig::new(seqname);
            new_contig.add_record(record);
            self.add_contig(new_contig);
        }
    }
    fn process_staged_records(&mut self) {
        for contig in self.contigs.iter_mut() {
            contig.process_staged_records()
        }
    }
    fn find_transcript<'a>(&'a self, searcher: &mut GtfSearcher<'a>) -> Option<&Transcript> {
        self.contigs
            .iter()
            .find_map(|c| c.find_transcript(searcher))
    }
}
