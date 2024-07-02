import subprocess


def extract_sequence(fasta_file, chromosome, start, end):

    command = (
        "samtools faidx --fai-idx "
        f"{fasta_file}.fai {fasta_file} {chromosome}:{start}-{end}"
    )
    result = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        text=True
    )
    if result.returncode == 0:
        return result.stdout.split('\n', 1)[1].replace('\n', '')
    else:
        print("Error:", result.stderr)
        return None
