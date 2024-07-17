extern crate genefeatures;

#[cfg(test)]
mod tests {
    use std::path::Path;
    use genefeatures::make_gf_fasta;


    #[test]
    fn test_read_gene_file() {
        let filepath: &Path = Path::new("tests/data/wt_gf_make_input_test.csv");
        let genes: Vec<(Box<str>, Option<Box<str>>, Option<Box<str>>)> = make_gf_fasta::read_gene_file(filepath);
        assert_eq!(genes[0], (Box::from("ENSG00000252182"), None, None))
    }
}