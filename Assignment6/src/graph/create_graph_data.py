"""
Module for creating networkx graph objects.
"""

import sys
import time
import argparse
from typing import Any
from pathlib import Path
import multiprocessing as mp
import os
import pickle

import pandas as pd
import numpy as np

import pyspark
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StringType, ArrayType, DateType

import pyspark.sql.functions as F

root_dir = os.path.abspath(os.path.join(
                  os.path.dirname(__file__),
                  os.pardir))

sys.path.insert(0, root_dir)


from parallel_computing.server_side import ServerSide
from parallel_computing.client_side import ClientSide

__author__ = "Stijn Arends"
__version__ = "v.01"

class ArgumentParser:
    """
    Class to parse the input arguments.
    """

    def __init__(self) -> None:
        """Initializer"""
        self.parser = self._create_argument_parser()
        # Print help if no arguments are supplied and stop the program
        if len(sys.argv) == 1:
            self.parser.print_help(sys.stderr)
            sys.exit(1)
        self.arguments = self.parser.parse_args()
        self.check_arg_combi()

    @staticmethod
    def _create_argument_parser():
        """
        Create an argument parser.

        :returns
        --------
        parser - ArgumentParser
        """
        parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
            description="Python script that parses pubmed articles",
            epilog="Contact: stijnarend@live.nl")

        parser.version = __version__

        parser.add_argument("-d", '--data_dir', action="store",
                    dest="d", required=False, type=str,
                    help="Location of the parsed pubmed data in json format.")

        parser.add_argument("-n", action="store",
                           dest="n", required=False, type=int, default=1,
                           help="Number of peons per client.")

        parser.add_argument("-o", '--out-dir', action="store",
                    dest="o", required=False, type=str,
                    help="The output directory to store the results - required if server mode is selected")

        parser.add_argument("-p", '--port_number',dest='p',
                        help="The port number that will be used",
                        required=True, type=int)

        parser.add_argument("--host", dest="host",
            help="Hosts used, first input is set as the server host",
            required=True, nargs="?")

        parser.add_argument('-v',
            '--version',
            help='Displays the version number of the script and exitst',
            action='version')

        command_group = parser.add_mutually_exclusive_group(required=True)
        command_group.add_argument('-s', dest='s',
            help='Server mode, can\'t be used together with -c',
            action='store_true')

        command_group.add_argument('-c', dest='c',
            help='Client Mode, can\'t be used together with -s',
            action='store_true')

        return parser

    def get_argument(self, argument_key: str) -> Any:
        """
        Method to get an input argument.

        :parameters
        -----------
        argument_key - str
            Name of command line argument.

        :returns
        --------
        value - Any
            Value of a command line argument
        """
        if self.arguments is not None and argument_key in self.arguments:
            value = getattr(self.arguments, argument_key)
        else:
            value = None
        return value

    def check_arg_combi(self) -> None:
        """
        Check if the data and output directory are specified when also
        server mode is specified.
        """
        if self.get_argument('s') and not self.get_argument('d'):
            self.parser.print_help(sys.stderr)
            sys.exit("\nIf server mode has been selected please also provide the "\
                "path to the data using the -d, or --data-dir flags.")

        if self.get_argument('s') and not self.get_argument('o'):
            self.parser.print_help(sys.stderr)
            sys.exit("\nIf server mode has been selected please also provide the "\
                "path to the output folder using the -o, or --output-dir flags.")


def create_pyspark_df(file_path):
    """
    Create a pyspark dataframe from multiple json files.

    :parameters
    -----------
    file_path - str
        Directory containg json files used to create dataframe
    
    :returns
    --------
    df - pyspark DataFrame
        A pyspark dataframe
    sc - SparkContext
    spark - SparkSession
    """
    schema = StructType() \
        .add('pmid', StringType()) \
        .add('language', StringType()) \
        .add('author', StringType()) \
        .add('title', StringType()) \
        .add('co_authors', ArrayType(StringType())) \
        .add('journal', StringType()) \
        .add('key_words', ArrayType(StringType())) \
        .add('doi', StringType()) \
        .add('pmc', StringType()) \
        .add('publish_date', DateType()) \
        .add('ref_ids', ArrayType(StringType())) \
        .add("ref_authors", ArrayType(ArrayType(StringType()))) \
        .add('ref_titles', ArrayType(StringType())) \
        .add('ref_type', StringType())

    conf = pyspark.SparkConf().setAll([('spark.executor.memory', '128g'),
                                ('spark.master', 'local[16]'),
                                ('spark.driver.memory', '128g')])
    sc = SparkContext(conf=conf)
    sc.getConf().getAll()
    spark = SparkSession(sc)
    df = spark.read.option("multiline","true").schema(schema) \
        .json(file_path)
    return df, sc, spark

def create_adj_list(data, max_refs, out_adjlist_file) -> None:
    """
    Create an adjacency list using pandas dataframe.
    Where the first column contains the source and the rest
    of the columns contain the targets.

    :parameters
    -----------
    data - pd.DataFrame
        Data frame containing parsed pubmed articles
    max_refs - int
        The maximum number of references out of all the articles.
    out_adjlist_file - Path
        Location of output file
    """
    adjlist = pd.DataFrame(data['ref_ids'].values.tolist())
    adjlist.insert(0, 'pmid', data['pmid'].values)

    n_cols = len(adjlist.columns) 
    if not n_cols == max_refs + 1:
        to_add = (max_refs + 1) - n_cols
        for i in range(to_add):
            extra_col = pd.Series([np.NaN for _ in range(len(adjlist))])
            adjlist[f'{n_cols + i + 1}'] = extra_col

    adjlist = adjlist.astype('Int64')

    adjlist.to_csv(out_adjlist_file, sep=' ', header=False, index=False, mode='a')


def save_no_ref_nodes(data, out_node_list) -> None:
    """
    Append PMIDs(nodes) to a text file.

    :parameters
    -----------
    data - pd.DataFrame
        Data frame containing parsed pubmed articles
    out_node_list - Path
        Location of output file
    """
    nodes = data.pmid.values

    with open(out_node_list, 'ab') as fh:
        np.savetxt(fh, nodes, fmt='%4.0f', newline=' ')


def get_attribute_info(df) -> dict:
    """
    Get the attributes info for each article(node) by
    setting the PMID column as index and transforming the
    resulting data frame to a dictionary.

    :parameters
    -----------
    df - pd.DataFrame
        Parsed pubmed articles

    :returns
    --------
    attributes - dict
        The article(node) attributes
    """
    attr_df = df.copy()
    attr_df.drop_duplicates(subset=["pmid"], keep='last', inplace=True)
    attributes = attr_df.set_index('pmid').to_dict('index')
    return attributes


def read_pickl(file) -> Any:
    """
    Read in a pickle file.
    """
    with open(file, 'rb') as fh:
        data = pickle.load(fh)

    return data

def write_to_pickl(data, file) -> None:
    """
    Write a python object to a pickle file.
    """
    with open(file, 'wb') as fh:
        pickle.dump(data, fh)

def parse_graph_data(file, max_refs, out_adjlist_file, out_node_list, out_attr_dir) -> None:
    """
    Parse the article data to create data that can be used to create
    graphs; adjenceny list, attributes info, and list of nodes with no references(
    or edges).

    """
    df = pd.read_json(file)

    ref_ids = df[df.ref_type == "pmid"]
    if not len(ref_ids) == 0:
        data = ref_ids[["pmid", "ref_ids"]]
        create_adj_list(data, max_refs, out_adjlist_file)

    no_refs = df[df.ref_type != "pmid"]
 
    if not len(no_refs) == 0:
        save_no_ref_nodes(no_refs, out_node_list)

    attributes_data = get_attribute_info(df)
    write_to_pickl(attributes_data, out_attr_dir / (file.stem + ".pkl"))


def make_output_dir(path: Path) -> None:
    """
    Create a directory (if it does not exsit yet) to store the
    data.

    :parameter
    ----------
    path - Path
        Path to a directory

    :Excepts
    --------
    FileExistsError
        The directory already exists
    """
    try:
        path.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(f"[{make_output_dir.__name__}] {path} already exists.")


def combine_attribute_files(attributes_path):
    files = list(attributes_path.glob("*.pkl"))

    all_attributes = {}
    # Parse all files:
    for file in files:
        print(f"Parsing: {file}")
        data = read_pickl(file)
        all_attributes.update(data)
        os.remove(file)
    return all_attributes


def main():
    """
    main
    
    If this does not finish in time, try again and write out file names that have been processed
    to a log file and run it again exclusing all those files.
    """

    start_time = time.time()
    cla_parser = ArgumentParser()
    # cla_validator = CLIArgValidator()
    data_dir = cla_parser.get_argument('d')

    if data_dir:
        json_files = data_dir + "*.json"
        data_dir = Path(data_dir)
        data_files = list(data_dir.glob("*.json"))
        # cla_validator.validate_input_file(data_dir)

        print("Reading data:\n")
        spark_df, sc, spark = create_pyspark_df(json_files)
        print(spark_df.printSchema())

        print("Filter data:\n")
        pmid_refs = spark_df.filter(spark_df.ref_type == "pmid")
        pmid_refs.show()

        full_df = pmid_refs.select("pmid", 'ref_ids')

        print("Calculate max size:")
        max_size_full = full_df.select(F.size('ref_ids')).agg({'size(ref_ids)': 'max'}).take(1)[0][0]
        print(f"Max size: {max_size_full}")

        sc.stop()

    output_dir = cla_parser.get_argument('o')
    if output_dir:
        output_dir = Path(output_dir)
        output_attr = output_dir / "attribute_data/partitions_data/"
        make_output_dir(output_dir)
        make_output_dir(output_attr)
        out_adjlist_file = output_dir / "adjlist_articles.csv"
        out_node_list = output_dir / "nodes_no_refs.txt"
        out_attr_file = output_dir / "attribute_data/all_attributes.pkl"

    n_peons = cla_parser.get_argument('n')
    port = cla_parser.get_argument('p')
    host = cla_parser.get_argument('host')

    server_mode = cla_parser.get_argument('s')
    client_mode = cla_parser.get_argument('c')

    if server_mode:
        server_side = ServerSide(ip_adress=host, port=port,
            auth_key=b'whathasitgotinitspocketsesss?',
            poison_pill="MEMENTOMORI")
        server = mp.Process(target=server_side.run_server, args=(parse_graph_data,
            data_files, max_size_full, out_adjlist_file, out_node_list,
            output_attr)) # [:2]
        server.start()
        time.sleep(1)
        server.join()

    if client_mode:
        print('Selected client mode.')
        client_side = ClientSide(ip_adress=host, port=port,
            auth_key=b'whathasitgotinitspocketsesss?',
            poison_pill="MEMENTOMORI")
        client = mp.Process(target=client_side.run_client, args=(n_peons,))
        client.start()
        client.join()

    if not client_mode:
        all_attributes = combine_attribute_files(output_attr)
        print(len(all_attributes))
        write_to_pickl(all_attributes, out_attr_file)

    elapsed_time = time.time() - start_time
    days = 0
    if elapsed_time >= 86400:
        days = int(elapsed_time / 86400)
    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print(f"\n---- Time: {days}:{elapsed} ----")


# python create_adjaceny_list.py -d /commons/dsls/dsph/2022/final_parsed_articles/ --host assemblix2012 -p 4235 -s -o /commons/dsls/dsph/2022/graph_data/

if __name__ == "__main__":
    main()
