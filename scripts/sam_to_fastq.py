"""Convert samfile to fastq for mrsfast mapping."""

from __future__ import print_function
from __future__ import division

import argparse

def write_fastq(seq, outfile):
    name = "@0"
    qual = "?" * len(seq)
    print(name, seq, "+", qual, sep="\n", file=outfile)

def read_pass(flag, exclude_any):
    """Return whether read passes based on filters in exclude_any.
       Flag and exclude_any must be integers."""
    return flag & exclude_any == 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--min_length", type=int, default=36, help="Minimum length of reads to include (Default: %(default)s)")
    parser.add_argument("--offset", type=int, default=0, help="Bases to offset start of read (Default: %(default)s)")
    parser.add_argument("--exclude_any", type=int, default=3840, help="Exclude reads with specified flags by bitwise AND")

    args = parser.parse_args()

    with open(args.input, "r") as input, \
         open(args.output, "w") as output:
        for line in input:
            if line.startswith("@"):
                continue

            entries = line.rstrip().split()

            seq = entries[9]

            seq = seq[args.offset:]
            if len(seq) >= args.min_length and (int(entries[1]) & args.exclude_any == 0):
                write_fastq(seq, output)
