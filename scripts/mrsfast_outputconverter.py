#!/usr/bin/env python3
# Converts the output of mrsfast into a simple position and edit distance file for read counting

import sys
import gzip
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")
    parser.add_argument("--compress", action="store_true", help="Compress output with gzip")

    args = parser.parse_args()

    if args.compress:
        func = gzip.GzipFile
        mode = "wb"
    else:
        func = open
        mode = "w"

    with open(args.infile, "r") as infile, \
         func(args.outfile, mode) as outfile:

        for line in infile:
            #Skip sam header and mrsfast logging
            if line.startswith("@") or line.startswith("|") or line.startswith("-"):
                continue
            var = line.rstrip().split("\t")
            try:
                #checks for correctly formatted line with edit distance in form of NM:i:0/1/2
                editdistance = var[11].split(':')
                if editdistance[0].startswith("NM"):
                    outstring = "{}\t{}\t{}\n".format(var[2], int(var[3]), int(editdistance[2]))
                if args.compress:
                    outstring = outstring.encode('ascii')
                outfile.write(outstring)
            except IndexError as e:
                if line.startswith("ERR:"):
                    print("Mrsfast", line, file=sys.stderr)
                    sys.exit(1)
                if line.startswith("Total"):
                    print(line, file=sys.stderr)
                else:
                    print("IndexError:", e, line, file=sys.stderr)
