use std::ops::{Div, Range};
use std::str::FromStr;
use std::collections::HashMap;
use crate::gtf_record::GtfRecord;
use crate::gtf_tree::Transcript;


#[derive(Debug, PartialEq)]
pub enum Strand {
    Pos,
    Neg,
    NotSpecified
}
impl FromStr for Strand {
    type Err = ();
    fn from_str(input: &str) -> Result<Strand, Self::Err> {
        match input {
            "+" => Ok(Strand::Pos),
            "-" => Ok(Strand::Neg),
            "." => Ok(Strand::NotSpecified),
            _ => Ok(Strand::NotSpecified)
        }
    }
}


///// HELPER FUNCTIONS FOR SEQUENCE INDEX GENERATION  //////

// Gathers variables needed to generate sequence index
fn get_index_variables(transcript: &mut Transcript) -> (u64, u64, Strand, &mut Vec<GtfRecord>) {

    let strand: Strand = Strand::from_str(transcript.record.strand.as_ref())
        .expect("Cannot generate strand enum from transcript");
    let cds: &mut Vec<GtfRecord> = &mut transcript.cds_records.records;

    // find start and end of entire transcript accounting for strand 
    let (start, end): (u64, u64) = match strand {
        Strand::Pos | Strand::NotSpecified => {
            cds.sort_by_key(|rec| rec.start);
            (transcript.record.start, transcript.record.end)
        },
        Strand::Neg => {
            cds.sort_by_key(|rec| rec.end);
            (transcript.record.end, transcript.record.start)
            }
        };
    (start, end, strand, cds)
} 


// Gets start and end values from a CDS record, accounts for strand
fn get_interval_indices(strand: &Strand, rec: &GtfRecord) -> (u64, u64) {
    match strand {
        Strand::Pos | Strand::NotSpecified => (rec.start, rec.end + 1),
        Strand::Neg => (rec.end, rec.start - 1)
    }
}

// Adds elements with correct prefix and index to the mutation index
fn extend_mutation_index(
    mut_idx: &mut Vec<String>, start: u64, end: u64, prefix: String
) {
    for i in start..end {
        mut_idx.push(format!("{}{}", prefix, i))
    }
}


// Gerates index values for pre-start UTR and the first CDS
fn process_pre_tss_utr_and_first_cds(
    mutation_index: &mut Vec<String>,
    strand: &Strand,
    start: u64,
    mut cds_idx: u64,
    rec: Option<&mut GtfRecord>
) -> (u64, u64) {

    let (i_start, i_end): (u64, u64) = get_interval_indices(
        strand, rec.expect("No record found!")
    );
    // Process pre tss utr index
    extend_mutation_index(
        mutation_index,
        start.abs_diff(i_start),
        0, 
        "-".to_string()
    );
    // process first cds
    let cds_len =  i_end.abs_diff(i_start);
    extend_mutation_index(
        mutation_index,
        1,
        cds_len,
        "".to_string() 
    );
    cds_idx += cds_len;
    (i_end, cds_idx)
}


// Intronic portions of the index are relative to the nearest CDS, so they are processed in halves
fn process_intron_halves(mutation_index: &mut Vec<String>, i_start: u64, last_cds_end: u64, cds_idx: u64) {

    let intron_len = i_start
        .abs_diff(last_cds_end)
        .checked_add(1)
        .expect("Integer overflow occured!");
    let remainder: u64 = intron_len % 2;
    let mid_point: u64 = (intron_len + 1).div(2); 
    let cds_ref: u64 = cds_idx.abs_diff(1);

    // process first half of intron
    extend_mutation_index(
        mutation_index,
        1,
        mid_point.abs_diff(remainder),
        format!("{cds_ref}+")
    );
    //process second half of intron 
    extend_mutation_index(
        mutation_index,
        mid_point,
        0,
        format!("{cds_idx}-")
    );
}


// Iterate through CDS records generating index values as it goes, checks for and processes introns
fn process_cds_and_introns<'a>(
    mutation_index: &mut Vec<String>,
    cds: &mut impl Iterator<Item=&'a mut GtfRecord>,
    mut cds_idx: u64,
    strand: &Strand,
    mut last_cds_end: u64,
) {
    for rec in cds {

        let (i_start, i_end): (u64, u64) = get_interval_indices(strand, rec);
        // Check For intron
        let intron_found: bool = match strand {
            Strand::Pos | Strand::NotSpecified => i_start.saturating_sub((-1 as i64) as u64) != last_cds_end,
            Strand::Neg => i_start.saturating_sub(1) != last_cds_end
        };
        if intron_found {
            process_intron_halves(mutation_index, i_start, last_cds_end, cds_idx)
        }
        // process cds
        let cds_len = i_end.abs_diff(i_start); 
        extend_mutation_index(
            mutation_index,
            cds_idx,
            cds_idx.checked_add(cds_len)
                .expect("integer overflow occured"),
            "".to_string() 
        );
        last_cds_end = i_end;   
        cds_idx += cds_len;
    }
}


//// SEQUENCE INDEX DEFINITION AND IMPLEMENTATION ////


#[derive(Debug)]
pub struct SeqIdx {
    pub genomic_index: HashMap<u64, u64>,
    pub mutation_index: HashMap<String, u64>,
    pub mutation_log: Vec<(String, u64, u64, u64)>
}
impl SeqIdx {

    pub fn new(transcript: &mut Transcript) -> Self {

        let (start, end, strand, cds): (u64, u64, Strand, &mut Vec<GtfRecord>) = get_index_variables(transcript);
        let genomic_index: HashMap<u64, u64> = Self::init_genomic_index(start, end);
        let mutation_index: HashMap<String, u64> = Self::init_mutation_index(start, end, &strand, cds);

        Self {
            genomic_index,
            mutation_index,
            mutation_log: Vec::new(),
        }
    }

    fn init_genomic_index(start: u64, end: u64) -> HashMap<u64, u64> {
            Range {start: start, end: end}
                .into_iter()
                .zip(Range {start: 0, end: end - start} .into_iter())
                .collect()
    }

    fn init_mutation_index(start: u64, end: u64, strand: &Strand, cds: &mut Vec<GtfRecord>) -> HashMap<String, u64> {

        // init empty index and prepare intervals for iteration
        let mut mutation_index: Vec<String> = Vec::with_capacity(start.abs_diff(end) as usize);
        let mut cds: std::slice::IterMut<GtfRecord> = cds.into_iter();

        //Process pre-start utr and first cds
        let (last_cds_end, cds_idx): (u64, u64) = process_pre_tss_utr_and_first_cds(
            &mut mutation_index,
            strand,
            start,
            1,
            cds.next()
        );
        //process introns and remaining cds intervals
        process_cds_and_introns(
            &mut mutation_index,
            &mut cds,
            cds_idx,
            strand,
            last_cds_end
        );
        // process post stop UTR
        extend_mutation_index(
            &mut mutation_index,
            1,
            end.abs_diff(last_cds_end),
            "*".to_string() 
        );
        //return mutation index
        mutation_index.into_iter()
            .zip( Range { start: 0, end: end.abs_diff(start) } .into_iter() )
            .collect()
    }
}
