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
    # with open(file, 'r') as f:
    groups = groupby(file, key=lambda x: not x.startswith(">"))
    # print(len(list(groups)))
    for k,v in groups:
        if k:
            count += 1
            #if count < 3:
            seq_list = list(v)
            # print(f"Value: {list(v)}")
            # seq = "".join(map(str.rstrip, list(v),""))
            seq = "".join(map(str.rstrip, seq_list))
            seqs.append(len(seq))

    return sorted(seqs)[::-1]


def calculate_N50(list_of_lengths):
    """Calculate N50 for a sequence of numbers.
 
    Args:
        list_of_lengths (list): List of numbers.
 
    Returns:
        float: N50 value.
 
    """
    tmp = []

    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
    # print(tmp)
 
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
 
    return median


# def test(line_lengths):
#     print(len(line_lengths))
#     print(len(line_lengths) % 2 == 0)

#     print(f"Count: {line_lengths.count(2050)}")
#     counts = {}

#     for line in set(line_lengths):
#         counts[line] = line_lengths.count(line) / len(line_lengths)

#     # print(counts)
#     dic2=dict(sorted(counts.items(),key= lambda x:x[1])[::-1])
    
#     print(dic2)

#     print(sum(dic2.values()))

#     # val = 0

#     # while val < 0.5:

#     #     for n_nucs, weight in dic2.items():
#     #         val += weight
#     odd = len(line_lengths) % 2 != 0

#     # calc_w_medium(dic2, odd)


def calc_N50(line_seqs):
    """
    Derived from : https://www.molecularecologist.com/2017/03/29/whats-n50/
    """
    # print(line_seqs)

    # total_length = sum(line_seqs)
    # print(f"Total assembly length: {total_length}")

    # val = 0
    # for size in line_seqs:
    #     val += size
    #     if val >= total_length /2:
    #         print(f"N50: {size}")
    #         break


    cs = np.cumsum(line_seqs)
    half_point = sum(line_seqs) / 2
    index = np.searchsorted(cs, half_point)
    n50 = line_seqs[index]
    # print(f"N50: {n50}")
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
    # contigs = sys.argv[1]
    # seq_lengths = read_file(contigs)
    file = sys.stdin
    # print(len(file))
    seq_lengths = read_file(file)


    n50 = calculate_N50(seq_lengths)
    # print(f"N50: {n50}")

    # test(seq_lengths)

    n50 = calc_N50(seq_lengths)

    # sys.stdout.write(f"{n50}\n")
    sys.stdout.write(f"{kmr_size},{n50}\n")

    # sys.stdout.write(f"blablbabla\n")


if __name__ == "__main__":
    main() # cat /commons/dsls/dsph/2022/velvet_out/contigs.fa | python3 assignment4.py 