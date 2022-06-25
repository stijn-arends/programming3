#!/usr/bin/env python3
from __future__ import division
import sys
from itertools import groupby


import numpy as np
import argparse
import os


def read_file(file):
    """
    https://stackoverflow.com/questions/29805642/learning-to-parse-a-fasta-file-with-python
    """
    seqs = []
    count = 0
    groups = groupby(file, key=lambda x: not x.startswith(">"))
    for k,v in groups:
        if k:
            count += 1
            seq_list = list(v)
            seq = "".join(map(str.rstrip, seq_list))
            seqs.append(len(seq))

    return sorted(seqs)[::-1]

def calc_N50(line_seqs):
    """
    Derived from : https://www.molecularecologist.com/2017/03/29/whats-n50/
    """
    cs = np.cumsum(line_seqs)
    half_point = sum(line_seqs) / 2
    index = np.searchsorted(cs, half_point)
    n50 = line_seqs[index]
    return n50


def get_arguments():
    """
    Create an argument parser.

    :returns
    --------
    parser - ArgumentParser
    """
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
        description="Plot some data",
        epilog="Contact: stijnarend@live.nl")

    parser.add_argument("-kmr", "--kmr_size", action="store", type=int, nargs=1,
            help="KMR size used for assembly", dest="kmr", required=True)

    parser.add_argument("-o", "--output", dest='out',
        action="store", type=str,
            help="Name of output file")

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    kmr_size = args.kmr[0]
    file = sys.stdin
    seq_lengths = read_file(file)

    n50 = calc_N50(seq_lengths)

    sys.stdout.write(f"{kmr_size},{n50}\n")



if __name__ == "__main__":
    main() # cat /commons/dsls/dsph/2022/velvet_out/contigs.fa | python3 assignment4.py 