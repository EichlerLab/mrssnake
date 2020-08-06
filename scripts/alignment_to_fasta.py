#!/bin/env python3

import pysam
import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument("--input", "-i", type=str, required=True, help="Input CRAM file")
parser.add_argument("--chunk", "-c", type=int, required=False, default=36, help="Size of reads to chunk sequence into")
parser.add_argument("--output", "-o", type=str, required=True, help="Output FASTA file")
parser.add_argument("--compression", "-p", type=str, required=True, help="Compression of file rb for bam, rc for cram")

args = parser.parse_args()

samfile = pysam.AlignmentFile(args.input, mode=args.compression, check_sq=False)
with open(args.output, 'w') as outFile:
	for line in samfile.fetch():
		if ("H" not in str(line.cigarstring)) and ((line.flag & 0xc00) == 0):
			sequence = line.query_sequence
			for i in range(0, len(sequence), args.chunk):
				if i+args.chunk < len(sequence):
					outFile.write(">0\n"+sequence[i:i+args.chunk]+'\n')
		else:
			continue

