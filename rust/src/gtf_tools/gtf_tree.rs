use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use gtf_record::GtfRecord;
use gtf_searcher::GtfSearcher;


#[derive(Debug)]
pub enum NodeType {
    Contig(Contig),
    Gene(Gene),
    Transcript(Transcript),
    Cds(Cds),
    NonCds(NonCds),
}

pub trait Node {
    fn add_record(&mut self, record: GtfRecord) {}

}

#[derive(Debug)]
pub struct Cds {
    pub records: Vec<GtfRecord>,
}
impl Cds {
    fn new () -> Self {
        Cds {
            records: Vec::new()
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
            records: Vec::new()
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
    pub non_cds_records: NonCds
}

impl Transcript {

    fn new(record: GtfRecord) -> Self {
        Transcript {
            record,
            cds_records: Cds::new(),
            non_cds_records:  NonCds::new()
        }
    }
}
impl Node for Transcript {
    fn add_record(&mut self, record: GtfRecord) {
        match record.feature.as_ref() {
            "CDS" => self.cds_records.add_record(record),
            _ => self.non_cds_records.add_record(record)
        }
    }
}

#[derive(Debug)]
pub struct Gene {
    pub record: GtfRecord,
    pub transcripts: Vec<Transcript>,
}
impl Gene {

    fn new(record: GtfRecord) -> Self {
        Gene {
            record,
            transcripts: Vec::new()
        }
    }

    fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts.push(transcript)
    }
}
impl Node for Gene {

    fn add_record(&mut self, record: GtfRecord) {
        match record.feature.as_ref() {
            "transcript" => self.add_transcript(Transcript::new(record)),
            _ => for t in self.transcripts.iter_mut() {
                if record.transcript_id == t.record.transcript_id {
                    t.add_record(record);
                    break
                }
            }
        }
    }
}

#[derive(Debug)]
pub struct Contig {
    pub genes: Vec<Gene>,
    pub name: String
}
impl Contig {

    fn new(name: String) -> Self {
        Contig {name, genes: Vec::new()}
    }

    fn add_gene(&mut self, gene: Gene) {
        self.genes.push(gene)
    }
}
impl Node for Contig {
    fn add_record(&mut self, record: GtfRecord) { 
        match record.feature.as_ref() {
            "gene" => self.add_gene(Gene::new(record)),
            _ => for g in self.genes.iter_mut() {
                if g.record.gene_id == record.gene_id {
                    g.add_record(record);
                    break
                }
            }
        }
    }
}


#[derive(Debug)]
pub struct GtfTree {
    pub contigs: Vec<Contig>,
    pub genome_build: Option<String>,
    pub genome_version: Option<String>,
    pub genome_date: Option<String>,
    pub genome_build_accension: Option<String>,
    pub genebuild_last_updated: Option<String>
}
impl GtfTree {

    fn new() -> Self {
        GtfTree {
            contigs: Vec::new(),
            genome_build: None,
            genome_version: None,
            genome_date: None,
            genome_build_accension: None,
            genebuild_last_updated: None
        }
    }

    fn add_contig(&mut self, contig: Contig) {
        self.contigs.push(contig)
    }

    fn add_metadata(&mut self, key: &str, val: &str) {
        match key {
            "genome-build" => self.genome_build = Some(val.to_string()),
            "genome-version" => self.genome_version = Some(val.to_string()),
            "genome-date" => self.genome_date = Some(val.to_string()),
            "genome-build-accession" => self.genome_build_accension = Some(val.to_string()),
            "genebuild-last-updated" => self.genebuild_last_updated = Some(val.to_string()),
            _ => {}
        }
    }

    fn parse_gtf_file(gtf_path: &Path) -> Self {

        let gtf: File = match File::open(gtf_path) {
            Ok(gtf) => gtf,
            Err(gtf) => panic!(gtf)
        };
        let mut gtf_tree: GtfTree = GtfTree::new();

        for line in io::BufReader::new(gtf).lines() {
            let mut to_record = match line {
                Ok(l) => l,  
                Err(_) => continue
            };
            if to_record.starts_with("#!") {
                let mut iter: std::str::Split<&str> = to_record
                    .strip_prefix("#!")
                    .expect("invalid metadata line!")
                    .split(" ");
                let key: &str = iter.next().expect("invalid metadata key");
                let val: &str = iter.next().expect("invalid metadata value");
                gtf_tree.add_metadata(key, val);
            } else {
                let record = GtfRecord::from_gtf_line(&to_record);
            }
        }
        gtf_tree
    }
}
impl Node for GtfTree {
    fn add_record(&mut self, record: GtfRecord) {
    }
}