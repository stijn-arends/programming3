"""
Module for mapping authors + title to their PMID.
"""


# IMPORTS
import glob
import ast
import time
import sys
import argparse
import os
from pathlib import Path
from typing import Any
import multiprocessing as mp
import json

import pandas as pd
import pyspark
from pyspark import SparkContext
from pyspark import SparkFiles
from pyspark.sql import SQLContext
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StringType, ArrayType, DateType

import pyspark.sql.functions as F
import pyspark.sql.types as T

from parallel_computing.server_side import ServerSide
from parallel_computing.client_side import ClientSide
# from server_side import ServerSide
# from client_side import ClientSide

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


# def create_pyspark_df(file_path):

#     schema = StructType() \
#         .add('pmid', StringType()) \
#         .add('language', StringType()) \
#         .add('author', StringType()) \
#         .add('title', StringType()) \
#         .add('co_authors', ArrayType(StringType())) \
#         .add('journal', StringType()) \
#         .add('key_words', ArrayType(StringType())) \
#         .add('doi', StringType()) \
#         .add('pmc', StringType()) \
#         .add('publish_date', DateType()) \
#         .add('ref_ids', ArrayType(StringType())) \
#         .add("ref_authors", ArrayType(ArrayType(StringType()))) \
#         .add('ref_titles', ArrayType(StringType())) \
#         .add('ref_type', StringType())

#     conf = pyspark.SparkConf().setAll([('spark.executor.memory', '128g'),
#                                 ('spark.master', 'local[16]'),
#                                 ('spark.driver.memory', '128g')])
#     sc = SparkContext(conf=conf)
#     sc.getConf().getAll()
#     # sc = SparkContext('local[16]')
#     spark = SparkSession(sc)
#     # df = spark.read.json(file)
#     df = spark.read.option("multiline","true").schema(schema) \
#         .json(file_path)
        
#     df = df.withColumn('row_id',F.monotonically_increasing_id())
#     return df, sc, spark

# def get_essential_data(spark_df):
#     essentials = spark_df.select('pmid', 'author', 'title')
#     # Remove dots from title
#     essentials = essentials.withColumn('title', F.regexp_replace('title', r"[\.]+", ""))

#     # Remove comma and whitespace beginning of title ^[, ]+
#     essentials = essentials.withColumn('title', F.regexp_replace('title', r"^[, ]+", ""))

#     # remove comma and whitespace at the end of title
#     essentials = essentials.withColumn('title', F.regexp_replace('title', r"[, ]+$", ""))

#     return essentials

def read_data_pandas(files):
    """
    It takes 3 minutes and 51 seconds to load.
    """
    df = pd.concat((pd.read_json(f) for f in files))
    return df


def map_author_to_pmid(file, essentials, output):
    file = Path(file)
    print(f"Processing file: {file.stem}")
    df = pd.read_json(file)
    # print(df)
    author_refs = df[df.ref_type == 'author']
    # print(author_refs)
    total = len(author_refs)
    count = 0
    for index, row in author_refs.iterrows():
        count+=1

        # print(f"\nProcessing sample: {count}")
        ref_authors = row['ref_authors']
        ref_titles = row['ref_titles']
        ref_titles = [title.replace(".", "").lstrip(",").lstrip().rstrip(",").rstrip() for title in ref_titles]
        pmid_list = []

        for _, (ref_title, ref_author) in enumerate(zip(ref_titles, ref_authors)):
            if ref_author:
                main_author = ref_author[0]
                exists = essentials[(essentials.author == main_author) & (essentials.title == ref_title)]
                # exists = essentials.filter((F.col('title') == ref_title) & (F.col('author') == main_author))
                pmid = ""
                if not exists.empty:
                    pmid = exists['pmid'].values[0]
                    
            else:
                pmid = ""
                
            pmid_list.append(pmid)
        df['ref_ids'].iloc[index] = pmid_list
        print(f"Parsed: {(count / total) * 100:.3f} %")

    json_file = output / (file.stem + "_mapped.json")

    json_data = json.loads(df.to_json(orient='records', date_format='iso'))

    with open(json_file, 'w') as f:
        f.write(json.dumps(json_data))


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

def main():
    """
    Main script.
    """
    start_time = time.time()
    cla_parser = ArgumentParser()
    # cla_validator = CLIArgValidator()
    data_dir = cla_parser.get_argument('d')

    if data_dir:
        data_dir = Path(data_dir)
        data_files = list(data_dir.glob("*.json"))
        # cla_validator.validate_input_file(data_dir)
        start_time_read_data = time.time()
        print("Reading the data.")
        df_pd = read_data_pandas(data_files)
        elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time_read_data))
        print(f"\n---- Time reading data: {elapsed} ----")

        print("\nPreparing the essential columns")
        essentials = df_pd.loc[:, ['pmid', 'author', 'title']]
        essentials["title"] = essentials["title"].str.replace('.', '', regex=False).str.lstrip(",")\
            .str.lstrip().str.rstrip(",").str.rstrip()

    output_dir = cla_parser.get_argument('o')
    if output_dir:
        output_dir = Path(output_dir)
        make_output_dir(output_dir)

    n_peons = cla_parser.get_argument('n')
    port = cla_parser.get_argument('p')
    host = cla_parser.get_argument('host')

    server_mode = cla_parser.get_argument('s')
    client_mode = cla_parser.get_argument('c')
    # data_dir = "/commons/dsls/dsph/2022/test/"

    # json_files = glob.glob(data_dir + "/*.json") # change this to Pathlib.glob()
    # print(json_files)
    # print(type(json_files))




    # essentials = pd.DataFrame({'pmid':['123456'], 'author': ['Stijn, A.'], 'title':'Some title'})
    # print(essentials)

    #     # ---- Test server client mode ------

    if server_mode:
        # out_dir = Path("/commons/dsls/dsph/2022/test")
        # out_dir.mkdir(exist_ok=True)
        server_side = ServerSide(ip_adress=host, port=port,
            auth_key=b'whathasitgotinitspocketsesss?',
            poison_pill="MEMENTOMORI")
        server = mp.Process(target=server_side.run_server, args=(map_author_to_pmid,
            data_files, essentials, output_dir)) # [:2]
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

    elapsed_time = time.time() - start_time
    days = 0
    if elapsed_time >= 86400:
        days = int(elapsed_time / 86400)
    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print(f"\n---- Time: {days}:{elapsed} ----")


if __name__ == "__main__":
    main()