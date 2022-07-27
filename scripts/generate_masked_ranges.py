#!bin/python

# Adapted from https://www.danielecook.com/generate-a-bedfile-of-masked-ranges-a-fasta-file/

import gzip
import sys

# This file will generate a BED file of the masked regions (Ns) of the FASTA file.

# STDIN or arguments
if len(sys.argv) > 1:
    # Check file type
    if sys.argv[1].endswith(".fa.gz") or sys.argv[1].endswith(".fna.gz") or sys.argv[1].endswith(".fasta.gz"):
        input_fasta = gzip.open(sys.argv[1], 'rt')
    elif sys.argv[1].endswith(".fa") or sys.argv[1].endswith(".fna") or sys.argv[1].endswith(".fasta"):
        input_fasta = open(sys.argv[1], 'rt')
    else:
        raise Exception("Unsupported File Type [supported '.fa', '.fna', '.fasta', '.fa.gz', '.fna.gz', '.fasta.gz']")
else:
    print("""
    \tUsage:\n\t\tgenerate_masked_ranges.py <fasta file | .fa or .fa.gz> <chrome find> <chrome replace>
    
    \t\t'Chrome find' and 'chrome replace' are used to find and replace the name of a chromosome. For example,
    \t\treplacing CHROMOSOME_I with chr1 can be accomplished by using the command as follows:
    \t\t\tpython generate_masked_ranges.py my_fasta.fa CHROMOSOME_ chr
    \t\tOutput is to stdout
    """)
    raise SystemExit

n, state = 0, 0  # line, character, state (0=Out of gap; 1=In Gap)
chrom, start, end = None, None, None

with input_fasta as f:
    for line in f:
        line = line.replace("\n", "")
        if line.startswith(">"):
            # Print end range
            if state == 1:
                print('\t'.join([chrom, str(start), str(n)]))
                start, end, state = 0, 0, 0
            n = 0  # Reset character
            chrom = line.split(" ")[0].replace(">", "")
            # If user specifies, replace chromosome as well
            if len(sys.argv) > 2:
                chrom = chrom.replace(sys.argv[2], sys.argv[3])
        else:
            for char in line:
                if state == 0 and char == "N":
                    state = 1
                    start = n
                elif state == 1 and char != "N":
                    state = 0
                    end = n
                    print('\t'.join([chrom, str(start), str(end)]))
                else:
                    pass

                n += 1  # First base is 0 in bed format.

# Print mask close if on the last chromosome.
if state == 1:
    print('\t'.join([chrom, str(start), str(n)]))
    start, end, state = 0, 0, 0
