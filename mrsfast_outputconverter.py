#!/usr/bin/env python3
# Converts the output of mrsfast into a simple position and edit distance file for read counting

import sys

myFile=sys.stdin

for line in myFile:
    if line.startswith("@"):
        continue
    line = line.rstrip()
    #skip sam header
    var = line.split('\t')
    try:
        #checks for correctly formatted line with edit distance in form of NM:i:0/1/2
        editdistance = var[11].split(':')
        if editdistance[0].startswith("NM"):
            print(var[2] + '\t' + str(int(var[3])) + '\t' + editdistance[2])
    except IndexError as e:
        print("IndexError:", e, file=sys.stderr)
