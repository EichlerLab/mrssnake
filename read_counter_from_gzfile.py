# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False
"""
Convert mrsfast sam output to matrix of read depth and read start depth values.
Uses sparse matrices to reduce memory footprint.
"""

from __future__ import print_function
from __future__ import division

import sys
import gzip
import argparse

from functools import total_ordering
import time

import re

import tables
import numpy as np
import pandas as pd
from tables import NoSuchNodeError

try:
    import create_depth_array
except ImportError as e:
    import pyximport
    pyximport.install()
    import create_depth_array

def create_array(contigs, contigs_file, max_edist):
    """Create global numpy array 'matrix' for given contig.
    shape = (contig length, nedists, 2)
    third dimension is for depth and starts, respectively.
    """

    contig_dat = pd.read_table(contigs_file, header=None, names=["contig", "length"])
    contig_dat.index = contig_dat.contig

    global matrix_dict
    matrix_dict = {}

    if contigs is None:
        contigs = contig_dat.contig
    for contig in contigs:
        contig_length = contig_dat.ix[contig_dat.contig == contig, "length"]
        matrix_dict[contig] = np.ndarray((int(contig_length), int(max_edist+1)), dtype=np.uint16)

def add_to_array(array, pos, edist, rlen=36):
    """Add read entry to contig. Assumes pos is 0-based.
    """

    end = pos + rlen
    array[pos, edist] += 1

def process_file(infile, contigs, max_edist, rlen=36, mode="tab"):
    contigs_string = "|".join(contigs)
    if mode == "tab":
        regex_full = re.compile("^({})\t([0-9]+)\t([0-9]+)".format(contigs_string))
    else:
        regex_full = re.compile("[^ @\t]+\t[0-9]+\t(%s)\t([0-9]+)\t.+NM:i:([0-9]+)" % contigs_string)
    nhits = {contig: 0 for contig in contigs}
    with gzip.open(infile, "rt") as handle:
        for line in handle:
            match = regex_full.match(line)
            if match is not None:
                contig, pos, edist = match.group(1,2,3)
                pos = int(pos) - 1 # Convert 1-based pos to 0-based
                edist = int(edist) 
                add_to_array(matrix_dict[contig], pos, edist, rlen=36)
                nhits[contig] += 1

    return nhits

def create_depth_contig(contig, read_len=36):
    """Takes contig name and uses global matrix_dict to get read start matrix.
    Uses that to create and return read depth matrix.
    """
    depth_count = np.zeros(shape=matrix_dict[contig].shape)

    # Matrix is contig_length x nedists

    for i in range(matrix_dict[contig].shape[0] - read_len):
        for j in range(matrix_dict[contig].shape[1]):
            if matrix_dict[contig][i, j] != 0:
                depth_count[i:i+read_len, j] += matrix_dict[contig][i, j]

    return depth_count

def write_to_h5(contig, depth_contig, fout_handle, chunksize=1000000):
    """Write counts (dictionary of contig matrices) to fout hdf5 file
    in increments of chunksize bases. Outfile is in wssd_out_file format.
    """
    try:
        group = fout_handle.get_node(fout_handle.root, "depthAndStarts_wssd")
    except NoSuchNodeError:
        group = fout_handle.create_group(fout_handle.root, "depthAndStarts_wssd")
    finally:
        contig_len, edists = matrix_dict[contig].shape
        carray_empty = tables.CArray(group,
                                     contig,
                                     tables.UInt16Atom(),
                                     (contig_len, edists, 2),
                                     filters=tables.Filters(complevel=1, complib="lzo")
                                    )

        nchunks = contig_len // chunksize
        if nchunks * chunksize < contig_len:
            nchunks += 1

        for i in range(nchunks):
            s = i * chunksize
            e = s + chunksize
            if e > contig_len:
                e = contig_len
            carray_empty[s:e, :, 0] = depth_contig[s:e, :]
            carray_empty[s:e, :, 1] = matrix_dict[contig][s:e, :]

        fout_handle.flush()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input file (sam or chr,pos,edist tab format)")
    parser.add_argument("outfile", help="HDF5 file with read depths and starts")
    parser.add_argument("contigs", nargs="+", help="Names of contigs")
    parser.add_argument("--contigs_file", required=True,
                        help="tab-delimited file with contig names and lengths")
    parser.add_argument("--mode", choices=["tab", "sam"], default="tab",
                        help="Input file format (Default: %(default)s, tab must also be gzipped)")
    parser.add_argument("--max_edist",
                        default=2,
                        type=int,
                        help="Maximum edit distance of input reads (default: %(default)s)"
                       )
    parser.add_argument("--read_length", default=36, type=int,
                        help="Length of input reads (default: %(default)s)")
    parser.add_argument("--log", default=sys.stderr, help="Path to log file. Default: sys.stderr")

    args = parser.parse_args()

    run_start = time.time()

    if args.log is not sys.stderr:
        logfile = open(args.log, "w")
    else:
        logfile = sys.stderr

    create_array(args.contigs, args.contigs_file, args.max_edist)

    nhits = process_file(args.infile, args.contigs, args.max_edist, args.mode)

    mp_end = time.time() - run_start

    print("Finished reading samfile in %d seconds." % mp_end, file=logfile)
    for contig in args.contigs:
        print("Contig %s had %d hits" % (contig, nhits[contig]), file=logfile)
    print("Writing to hdf5: %s" % args.outfile, file=logfile)
    
    with tables.open_file(args.outfile, mode="w") as h5file:
        for contig in args.contigs:
            depth_contig = np.zeros(shape=matrix_dict[contig].shape, dtype=np.uint16)
            create_depth_array.create_depth_array(matrix_dict[contig], depth_contig)
            write_to_h5(contig, depth_contig, h5file)

    run_end = time.time()
    total_runtime = run_end - run_start

    print("Counter: finished writing matrices to hdf5: %s" % args.outfile, file=logfile)
    print("Counter: total runtime: %d seconds" % total_runtime, file=logfile, flush=True)
    if logfile is not sys.stderr:
        logfile.close()
    sys.exit()
