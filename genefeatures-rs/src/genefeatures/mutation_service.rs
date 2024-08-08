use crate::seq_index::SeqIdx;
use bio_seq::prelude::*;
use lazy_static::lazy_static;
use regex::{Captures, Regex};

// Compile all regex strings statically
lazy_static! {
    static ref SUB: Regex = Regex::new(r"(\d+)([ACGT])>([ACGT])").unwrap();
    static ref PNTDEL: Regex = Regex::new(r"^(\d+)del([ACGT]*)").unwrap();
    static ref RNGDEL: Regex = Regex::new(r"^(\d+)_(\d+)del([ACGT]*)").unwrap();
    static ref INS: Regex = Regex::new(r"(\d+)_(\d+)ins([ACGT]+)").unwrap();
    static ref DUP: Regex = Regex::new(r"(\d+)_(\d+)dup([ACGT]*)").unwrap();
    static ref INV: Regex = Regex::new(r"(\d+)_(\d+)inv(\d*)").unwrap();
    static ref INDEL: Regex =
        Regex::new(r"^(\d+)_(\d+)delins([ACGT]+)|^(\d+)_(\d+)del([ACGT]+)ins([ACGT]+)").unwrap();
}

#[derive(Debug, Clone)]
pub enum MutRegex {
    Sub,
    PntDel,
    RngDel,
    Ins,
    Dup,
    Inv,
    InDel,
}
impl MutRegex {
    fn get_regex(&self) -> &'static Regex {
        match self {
            MutRegex::Sub => &SUB,
            MutRegex::PntDel => &PNTDEL,
            MutRegex::RngDel => &RNGDEL,
            MutRegex::Ins => &INS,
            MutRegex::Dup => &DUP,
            MutRegex::Inv => &INV,
            MutRegex::InDel => &INDEL,
        }
    }

    pub fn is_match(&self, mutation: &str) -> bool {
        self.get_regex().is_match(mutation)
    }
    pub fn captures<'a>(&self, mutation: &'a str) -> Option<Captures<'a>> {
        self.get_regex().captures(mutation)
    }
    pub fn get_recipe<'a>(&self, mutation: &'a str) -> Vec<&'a str> {
        if let Some(captures) = self.get_regex().captures(mutation) {
            captures
                .iter()
                .skip(1)
                .filter_map(|cap| cap.map(|s| s.as_str()))
                .collect()
        } else {
            vec![]
        }
    }
}

// do not change order of below, pntdel and rngdel regex also matches indels
// cannot do look around with rust regex
// match indel first as it will fail on pntdel and rngdel strs
const MUT_REGEX: [MutRegex; 7] = [
    MutRegex::InDel,
    MutRegex::Ins,
    MutRegex::PntDel,
    MutRegex::RngDel,
    MutRegex::Sub,
    MutRegex::Dup,
    MutRegex::Inv,
];

pub fn process_mutation(mutation: String, seq: Seq<Dna>, seq_idx: SeqIdx) -> Seq<Dna> {
    let mut_regex: Option<&MutRegex> = MUT_REGEX.iter().find(|mr| mr.is_match(&mutation));

    if let Some(mr) = mut_regex {
        match mr {
            MutRegex::InDel => dna_indel(seq, seq_idx, mr.get_recipe(&mutation)),
            MutRegex::Dup => dna_dup(seq, seq_idx, mr.get_recipe(&mutation)),
            MutRegex::Sub => dna_sub(seq, seq_idx, mr.get_recipe(&mutation)),
            MutRegex::PntDel => dna_point_del(seq, seq_idx, mr.get_recipe(&mutation)),
            MutRegex::RngDel => dna_range_del(seq, seq_idx, mr.get_recipe(&mutation)),
            MutRegex::Ins => dna_ins(seq, seq_idx, mr.get_recipe(&mutation)),
            MutRegex::Inv => dna_inv(seq, seq_idx, mr.get_recipe(&mutation)),
        }
    } else {
        seq
    }
}

fn mutate_sequence(
    seq: Seq<Dna>,
    start: usize,
    end: usize,
    seq_ref: Seq<Dna>,
    seq_alt: Seq<Dna>,
) -> Seq<Dna> {
    assert_eq!(
        &seq[start..end],
        &seq_ref,
        "Reference sequence does not match full sequence slice"
    );
    let mut mutated_seq: Seq<Dna> = seq[..start].to_owned();
    mutated_seq.extend(seq_alt.iter());
    mutated_seq.extend(seq[end..].to_owned().iter());
    mutated_seq
}

fn dna_sub(seq: Seq<Dna>, seq_idx: SeqIdx, recipe: Vec<&str>) -> Seq<Dna> {
    let start: usize = seq_idx.mutation_index[recipe[0]] as usize;
    let end: usize = start + 1;
    let seq_ref: Seq<Dna> = recipe[1].try_into().unwrap();
    let seq_alt: Seq<Dna> = recipe[2].try_into().unwrap();
    mutate_sequence(seq, start, end, seq_ref, seq_alt)
}

fn dna_point_del(seq: Seq<Dna>, seq_idx: SeqIdx, recipe: Vec<&str>) -> Seq<Dna> {
    println!("\nCALLING PNT DEL\n");
    let start: usize = seq_idx.mutation_index[recipe[0]] as usize;
    let end: usize = start + 1;
    let seq_ref: Seq<Dna> = match recipe[1] {
        "" => seq[start].to_owned(),
        _ => recipe[1].try_into().unwrap(),
    };
    let seq_alt: Seq<Dna> = "".try_into().unwrap();

    mutate_sequence(seq, start, end, seq_ref, seq_alt)
}

fn dna_range_del(seq: Seq<Dna>, seq_idx: SeqIdx, recipe: Vec<&str>) -> Seq<Dna> {
    let (start, end): (usize, usize) = (
        seq_idx.mutation_index[recipe[0]] as usize,
        (seq_idx.mutation_index[recipe[1]] + 1) as usize,
    );
    let seq_ref: Seq<Dna> = match recipe[2] {
        "" => seq[start..end].to_owned(),
        _ => recipe[2].try_into().unwrap(),
    };
    let seq_alt: Seq<Dna> = "".try_into().unwrap();

    mutate_sequence(seq, start, end, seq_ref, seq_alt)
}

fn dna_ins(seq: Seq<Dna>, seq_idx: SeqIdx, recipe: Vec<&str>) -> Seq<Dna> {
    let (start, end): (usize, usize) = (
        (seq_idx.mutation_index[recipe[0]] + 1) as usize,
        (seq_idx.mutation_index[recipe[1]]) as usize,
    );
    let seq_ref: Seq<Dna> = dna!("").into();
    let seq_alt: Seq<Dna> = recipe[2].try_into().unwrap();

    mutate_sequence(seq, start, end, seq_ref, seq_alt)
}

fn dna_indel(seq: Seq<Dna>, seq_idx: SeqIdx, recipe: Vec<&str>) -> Seq<Dna> {
    let (start, end): (usize, usize) = (
        seq_idx.mutation_index[recipe[0]] as usize,
        (seq_idx.mutation_index[recipe[1]] + 1) as usize,
    );
    let (seq_ref, seq_alt): (Seq<Dna>, Seq<Dna>) = if recipe.len() == 4 {
        (
            // seq ref is in recipe if length is 4
            recipe[2].try_into().unwrap(),
            recipe[3].try_into().unwrap(),
        )
    } else {
        (
            // no seq ref in receipt if length is == 3
            seq[start..end].try_into().unwrap(),
            recipe[2].try_into().unwrap(),
        )
    };
    mutate_sequence(seq, start, end, seq_ref, seq_alt)
}

fn dna_dup(seq: Seq<Dna>, seq_idx: SeqIdx, recipe: Vec<&str>) -> Seq<Dna> {
    let (start, end): (usize, usize) = (
        seq_idx.mutation_index[recipe[0]] as usize,
        (seq_idx.mutation_index[recipe[1]] + 1) as usize,
    );
    let seq_ref: Seq<Dna> = seq[start..end].to_owned();
    let mut seq_alt: Seq<Dna> = seq_ref.clone();
    seq_alt.extend(&seq_ref);

    mutate_sequence(seq, start, end, seq_ref, seq_alt)
}

fn dna_inv(seq: Seq<Dna>, seq_idx: SeqIdx, recipe: Vec<&str>) -> Seq<Dna> {
    let (start, end): (usize, usize) = (
        seq_idx.mutation_index[recipe[0]] as usize,
        (seq_idx.mutation_index[recipe[1]] + 1) as usize,
    );
    let seq_ref: Seq<Dna> = seq[start..end].to_owned();
    let seq_alt: Seq<Dna> = seq_ref.rev().collect();

    mutate_sequence(seq, start, end, seq_ref, seq_alt)
}
