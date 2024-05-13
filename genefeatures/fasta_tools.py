import os
import subprocess


def extract_seq(fasta, seqname, start, stop):

    cmd = f"samtools faidx {fasta} {seqname}:{start}-{stop}"
    result = subprocess.run(cmd, shell = True, capture_output = True, text = True)

    if result.returncode == 0:
        return result.stdout.split('\n', 1)[1].replace('\n', '')
    else:
        print("Error:", result.stderr)
        return None