#!/bin/env python3

import os
import argparse
import gzip
from functools import reduce
import pysam

# Definitions




## Intake arguments
parser = argparse.ArgumentParser()

parser.add_argument("--file", "-f", type=str, required=True, help="Alignment file to be separated")
parser.add_argument("--index", "-i", type=str, required=True, help="Index along which the file will be split")
parser.add_argument("--compression", "-c", type=str, required=True, help="Compression type [rc: cram, rb : bam]")
parser.add_argument("--partitions", "-p", type=int, required=False, default=100, help="Number of partitions to split the file into")
parser.add_argument("--slice", "-s", type=int, required=False, help="Single partition number to extract")
parser.add_argument("--ref", "-r", type=str, required=True, help="Reference reads to which reads are aligned")
parser.add_argument("--output", "-o", type=str, required=False, default="/dev/stdout", help="Destination for parsed output")

args = parser.parse_args()

# Used to hold bin byte locations
coordList = []

# Reads in the byte location of bins of cram reads
with open(args.index) as f:
	possibleParts = [int(line.rstrip()) for line in f][0]

# Determines how many read bins can be contained per partition (rounds down)
perChunk = int(possibleParts/args.partitions)

# Used to keep track of how many bins per partition
binDirectory = [perChunk for totalBins in range(args.partitions)]

# Adjusts bin counts per partition so all bins are contained in partitions
binNum = 0
while sum(binDirectory) < possibleParts:
	binDirectory[binNum] += 1
	binNum += 1

i = args.slice-1

# Determines end partition
end = sum(binDirectory[0:i]) + binDirectory[i]
# Determines start partition
start = sum(binDirectory[0:i])

if args.compression == 'rc':
	samfile = pysam.AlignmentFile(args.file, mode=args.compression, reference_filename=args.ref)
else:
	samfile = pysam.AlignmentFile(args.file, mode=args.compression)

cram_chunk = pysam.AlignmentFile(args.output, "wb", template=samfile)

for i, read in enumerate(samfile.fetch()):
	if i < start:
		continue
	elif i > end:
		break
	else:
		cram_chunk.write(read)


cram_chunk.close()
samfile.close()





