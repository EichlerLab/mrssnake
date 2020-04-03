#!/bin/env python3

import os
import argparse
import gzip
from functools import reduce

# Definitions


### main function for copying byte-sized chunks to file
def readChunk(skip, size, mode):
	with open(args.cram, 'rb') as f:
		f.seek(skip)
		byte = f.read(size)
		with open(args.output, mode) as outfile:
			outfile.write(byte)

### Used on the last interval, to ensure proper EOF handling
def lastChunk(skip, mode):
	with open(args.cram, 'rb') as f:
		f.seek(skip)
		byte = f.read(1)
		with open(args.output, mode) as outfile:
			while byte:
				outfile.write(byte)
				byte = f.read(1)

## Intake arguments
parser = argparse.ArgumentParser()

parser.add_argument("--index", "-i", type=str, required=True, help="Index along which the file will be split")
parser.add_argument("--cram", "-c", type=str, required=True, help="CRAM to be separated")
parser.add_argument("--partitions", "-p", type=int, required=False, default=100, help="Number of partitions to split the file into")
parser.add_argument("--slice", "-s", type=int, required=False, help="Single partition number to extract")
parser.add_argument("--output", "-o", type=str, required=False, default="/dev/stdout", help="Destination for parsed output")

args = parser.parse_args()

# Used to hold bin byte locations
coordList = []

# Reads in the byte location of bins of cram reads
with open(args.index) as f:
	possibleParts = -1
	for line in f:
		line = line.rstrip()
		line = line.split(" ")
		possibleParts += 1
		coordList.append(int(line[2]))

# Last bin is empty EOF bin
EOFstart = coordList[-1]
# First bin starts after header
headerSize = coordList[0]
# Determines how many read bins can be contained per partition (rounds down)
perChunk = int(possibleParts/args.partitions)


binDirectory = []
# Used to keep track of how many bins per partition
for totalBins in range(args.partitions):
	binDirectory.append(perChunk)

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
# Determines where to start reading partitions
skip = coordList[start]

# Adjusts end if it continues past end of CRAM coordinates
if end > len(coordList)-1:
	readChunk(0, int(headerSize), 'wb')
	print("Finished with header")
	lastChunk(skip, 'ab')
	print("Finished with last bit")
else:
	# Calculates how many bytes to read in to file
	count = coordList[end]-coordList[start]
	readChunk(0, int(headerSize), 'wb')
	print("Finished with header")
	readChunk(int(skip), int(count), 'ab')
	print("Finished with body")
	lastChunk(int(EOFstart), 'ab')
	print("Finished with EOF")


