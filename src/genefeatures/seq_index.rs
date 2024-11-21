use crate::gtf_record::GtfRecordView;
use crate::gtf_tree::Transcript;
use std::cmp::Reverse;
use std::collections::HashMap;
use std::ops::Div;
use std::str::FromStr;

#[derive(Debug, PartialEq)]
pub enum Strand {
    Pos,
    Neg,
    NotSpecified,
}
impl FromStr for Strand {
    type Err = ();
    fn from_str(input: &str) -> Result<Strand, Self::Err> {
        match input {
            "+" => Ok(Strand::Pos),
            "-" => Ok(Strand::Neg),
            "." => Ok(Strand::NotSpecified),
            _ => Ok(Strand::NotSpecified),
        }
    }
}

///// HELPER FUNCTIONS FOR SEQUENCE INDEX GENERATION  //////

// Gathers variables needed to generate sequence index
fn get_index_variables(transcript: &Transcript) -> (u64, u64, Strand, Vec<GtfRecordView>) {
    let strand: Strand = Strand::from_str(transcript.record.strand.as_ref())
        .expect("Cannot generate strand enum from transcript");
    let mut cds_view: Vec<GtfRecordView> = transcript
        .cds_records
        .records
        .iter()
        .map(|rec| GtfRecordView::new(rec))
        .collect();

    // find start and end of entire transcript accounting for strand
    let (start, end): (u64, u64) = match strand {
        Strand::Pos | Strand::NotSpecified => {
            cds_view.sort_by_key(|rec| rec.start);
            (transcript.record.start, transcript.record.end + 1)
        }
        Strand::Neg => {
            cds_view.sort_by_key(|rec| Reverse(rec.end));
            (transcript.record.end, transcript.record.start - 1)
        }
    };
    (start, end, strand, cds_view)
}

// Gets start and end values from a CDS record, accounts for strand
fn get_interval_indices(strand: &Strand, rec: &GtfRecordView) -> (u64, u64) {
    match strand {
        Strand::Pos | Strand::NotSpecified => (rec.start.clone(), rec.end.clone() + 1),
        Strand::Neg => (rec.end.clone(), rec.start.clone() - 1),
    }
}

// Adds elements with correct prefix and index to the mutation index
fn extend_mutation_index(
    mut_idx: &mut Vec<String>,
    idx_iter: Box<dyn Iterator<Item = u64>>,
    prefix: String,
) {
    mut_idx.extend(idx_iter.map(|i| format!("{}{}", prefix, i)))
}

// Gerates index values for pre-start UTR and the first CDS
fn process_pre_tss_utr_and_first_cds(
    mutation_index: &mut Vec<String>,
    strand: &Strand,
    start: u64,
    mut cds_idx: u64,
    rec: Option<&GtfRecordView>,
) -> (u64, u64) {
    let (i_start, i_end): (u64, u64) =
        get_interval_indices(strand, rec.expect("not record view found"));
    // Process pre tss utr index
    let idx_iter: Box<dyn Iterator<Item = u64>> = Box::new((1..=start.abs_diff(i_start)).rev());
    extend_mutation_index(mutation_index, idx_iter, "-".to_string());
    // process first cds

    let cds_len = i_end.abs_diff(i_start);
    let idx_iter: Box<dyn Iterator<Item = u64>> = Box::new(1..=cds_len);
    extend_mutation_index(mutation_index, idx_iter, "".to_string());
    cds_idx += cds_len;
    (i_end, cds_idx)
}

// Intronic portions of the index are relative to the nearest CDS, so they are processed in halves
fn process_intron_halves(
    mutation_index: &mut Vec<String>,
    i_start: u64,
    last_cds_end: u64,
    cds_idx: u64,
) {
    let intron_len = i_start
        .abs_diff(last_cds_end)
        .checked_add(1)
        .expect("Integer overflow occured!");
    let remainder: u64 = intron_len % 2;
    let mid_point: u64 = (intron_len + 1).div(2);
    let cds_ref: u64 = cds_idx.abs_diff(1);

    // process first half of intron
    let iter_idx: Box<dyn Iterator<Item = u64>> = Box::new(1..mid_point.abs_diff(remainder));
    extend_mutation_index(mutation_index, iter_idx, format!("{cds_ref}+"));
    //process second half of intron
    let iter_idx: Box<dyn Iterator<Item = u64>> = Box::new((1..=mid_point).rev());
    extend_mutation_index(mutation_index, iter_idx, format!("{cds_idx}-"));
}

// Iterate through CDS records generating index values as it goes, checks for and processes introns
fn process_cds_and_introns<'a>(
    mutation_index: &mut Vec<String>,
    cds: core::slice::Iter<GtfRecordView>,
    mut cds_idx: u64,
    strand: &Strand,
    mut last_cds_end: u64,
) {
    for rec in cds {
        let (i_start, i_end): (u64, u64) = get_interval_indices(strand, rec);
        // Check For intron
        let intron_found: bool = match strand {
            Strand::Pos | Strand::NotSpecified => {
                i_start.saturating_sub((-1 as i64) as u64) != last_cds_end
            }
            Strand::Neg => i_start.saturating_sub(1) != last_cds_end,
        };
        if intron_found {
            process_intron_halves(mutation_index, i_start, last_cds_end, cds_idx)
        }
        // process cds
        let cds_len: u64 = i_end.abs_diff(i_start);
        let iter_idx: Box<dyn Iterator<Item = u64>> = Box::new(
            cds_idx
                ..cds_idx
                    .checked_add(cds_len)
                    .expect("iteger overflow occured"),
        );
        extend_mutation_index(mutation_index, iter_idx, "".to_string());
        last_cds_end = i_end;
        cds_idx += cds_len;
    }
}

//// SEQUENCE INDEX DEFINITION AND IMPLEMENTATION ////

#[derive(Debug)]
pub struct SeqIdx {
    pub genomic_index: HashMap<u64, u64>,
    pub sequence_index: HashMap<String, u64>,
    pub mutation_log: Vec<(String, u64, u64, u64)>,
}
impl SeqIdx {
    pub fn new(transcript: &Transcript) -> Self {
        let (start, end, strand, cds): (u64, u64, Strand, Vec<GtfRecordView>) =
            get_index_variables(transcript);
        let genomic_index: HashMap<u64, u64> = Self::init_genomic_index(start, end, &strand);
        let sequence_index: HashMap<String, u64> =
            Self::init_mutation_index(start, end, &strand, cds);

        Self {
            genomic_index,
            sequence_index,
            mutation_log: Vec::new(),
        }
    }

    fn init_genomic_index(start: u64, end: u64, strand: &Strand) -> HashMap<u64, u64> {
        let (genomic_pos, seq_len): (Box<dyn Iterator<Item = u64>>, u64) = match strand {
            Strand::Pos | Strand::NotSpecified => (Box::new(start..end), end.abs_diff(start)),
            Strand::Neg => (Box::new((end..=start).rev()), start.abs_diff(end)),
        };
        genomic_pos.zip(0..seq_len).into_iter().collect()
    }

    fn init_mutation_index(
        start: u64,
        end: u64,
        strand: &Strand,
        cds: Vec<GtfRecordView>,
    ) -> HashMap<String, u64> {
        // init empty index and prepare intervals for iteration
        let mut mutation_index: Vec<String> = Vec::with_capacity(start.abs_diff(end) as usize);
        let mut cds: std::slice::Iter<GtfRecordView> = cds.iter();

        //Process pre-start utr and first cds
        let (last_cds_end, cds_idx): (u64, u64) =
            process_pre_tss_utr_and_first_cds(&mut mutation_index, strand, start, 1, cds.next());
        //process introns and remaining cds intervals
        process_cds_and_introns(&mut mutation_index, cds, cds_idx, strand, last_cds_end);
        // process post stop UTR
        let iter_idx: Box<dyn Iterator<Item = u64>> = Box::new(1..end.abs_diff(last_cds_end));
        extend_mutation_index(&mut mutation_index, iter_idx, "*".to_string());
        //return mutation index
        mutation_index
            .into_iter()
            .zip(0..end.abs_diff(start))
            .into_iter()
            .collect()
    }
}
