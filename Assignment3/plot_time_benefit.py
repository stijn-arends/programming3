#!/usr/bin/env python3

"""
Simple script to plot some data
"""

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import os


def read_data(file: str) -> pd.DataFrame:
    """
    Read in some data.

    :parameters
    -----------
    file - str
        file name

    :returns
    --------
    df - pd.DataFrame
        df
    """
    df = pd.read_csv(file, sep=" ", names=['CPU', 'time'])
    return df


def plot_data(x, y, out_file) -> None:
    """
    Plot some data
    
    :parameters
    -----------
    x - array like
        X values
    y - array like
        Y values
    """
    plt.plot(x, y, 'o-')
    plt.xlabel('Number of threads used')
    plt.ylabel("Execution time (s)")
    plt.title("Comparing execution time for blastp using different number of CPUs")
    plt.savefig(out_file)


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

    parser.add_argument("in_file", action="store", type=str, nargs=1,
            help="input file")

    parser.add_argument("-o", "--output", dest='out',
        action="store", type=str,
            help="Name of output file")

    args = parser.parse_args()
    return args

def main():
    """main"""
    args = get_arguments()
    file = args.in_file[0]
    out_file = args.out

    # file = "timings.txt"
    df = read_data(file)

    x = df.CPU.values
    y = df.time.values

    plot_data(x, y, out_file)


if __name__ == "__main__":
    main()