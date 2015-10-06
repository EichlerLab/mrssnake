from __future__ import print_function
from __future__ import division

import sys
import argparse
import pysam
import numpy as np

SEND_BLOCK_SIZE = int(2e6)
READ_BLOCK_SIZE = int(10e6) # Chunk size for input

def chunk_reads(chr, start, end, bamfile, outfile, chunk_size, fifo):
    """This takes blocks of reads from a bam, then chunks them and sends them to mrsfast"""

    if chr == "unmapped":
        fetch_string = "*"
    elif start is None or end is None:
        fetch_string = chr
    else:
        fetch_string = "%s:%d-%d" % (chr, start + 1, end - 1)

    for l in pysam.view(bamfile, fetch_string):
        n_to_do = l.rlen // chunk_size
        for k in range(n_to_do):
            outfile.write(">0\n" + l.seq[k*chunk_size : k*chunk_size + chunk_size] + "\n")

    #Handle regions where there are no reads
    try:
        l
    except NameError:
        if args.fifo is not None:
            with open(fifo, "w") as fifo_handle:
                fifo_handle.write("ERROR: no reads for %s\n" % fetch_string)
        outfile.write("\n")
        sys.stderr.write("No reads for %s\n" % fetch_string)
    
def chunk_reads_old(chr, start, end, bamfile, outfile, chunk_size, dt):
    """This takes blocks of reads from a bam, then chunks them and sends them to mrsfast
       Code from psudmant's super_mapper.py"""
    reads_block = np.empty(SEND_BLOCK_SIZE, dtype = dt)

    curr_pos = 0
    for l in bamfile.fetch(chr, start, end):
        n_to_do = l.rlen // chunk_size
        for k in range(n_to_do):
            reads_block[curr_pos] = l.seq[k*chunk_size: k*chunk_size + chunk_size]
            curr_pos += 1
            if curr_pos == SEND_BLOCK_SIZE: 
                break
        if curr_pos == SEND_BLOCK_SIZE:
            curr_pos = 0
            outfile.write(reads_block)
            #READER_client_SEND_READS(reads_block,SEND_BLOCK_SIZE,proc_name,myrank,comm,icomm_to_RUNNER)
    
    if curr_pos!=0:
        small_reads_block = np.empty(curr_pos,dtype = dt)
        small_reads_block[:] = reads_block[:curr_pos]
        outfile.write(small_reads_block)
        #READER_client_SEND_READS(small_reads_block,curr_pos,proc_name,myrank,comm,icomm_to_RUNNER)

    ### make a new little block and copy stuff into it then send out
    #assert True

    print("FINISHED READING %s:%d-%d" % (chr,start,end))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bamfile")
    parser.add_argument("chr", help="Input bam contig to chunk. Use 'unmapped' to get all unmapped reads.")
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--outfile", default = "/dev/stdout")
    parser.add_argument("--chunk_size", default = 36, type = int)
    parser.add_argument("--fifo", help = "Path to fifo file")

    args = parser.parse_args()

    #bamfile = pysam.AlignmentFile(args.bamfile, "rb")
    outfile = open(args.outfile, "w")

    chunk_reads(args.chr, args.start, args.end, args.bamfile, outfile, args.chunk_size, args.fifo)

    #bamfile.close()
    outfile.close()
