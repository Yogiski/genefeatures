use std::process::{Command, Output};
use std::error::Error;


pub fn call_samtools_faidx(fasta: &str, chromosome: &str, start: u64, end: u64) -> Result<String, Box<dyn Error>> { 

    let faidx: String = format!("{fasta}.fai"); 
    let region: String = format!("{chromosome}:{start}-{end}");
    let args: [&str; 5]  = [
        "faidx", "--fai-idx", &faidx, fasta, &region
    ];
    let output: Output = Command::new("samtools")
        .args(args)
        .output()?;

    if output.status.success() {
        let res_string: String = String::from_utf8(output.stdout)?;
        let result: String = res_string.lines().skip(1).collect::<String>();
        Ok(result.trim().to_string())
    } else {
        // Command failed, get the error output
        let stderr: std::borrow::Cow<str> = String::from_utf8_lossy(&output.stderr);
        // Return the error output as a String wrapped in Err
        Err(format!("Command `samtools {} {} {} {}` failed with error: {}", 
                     args[0], args[2], args[3], args[4], stderr).into())
    }
}
