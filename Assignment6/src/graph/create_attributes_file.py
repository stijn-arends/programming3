#!/usr/bin/env python

"""Module"""

__author__ = 'Stijn Arends'
__version__ = 'v.01'

# IMPORTS
import sys
import os
import argparse
from pathlib import Path
from typing import Any

import pyspark
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StringType, ArrayType, DateType

from pyspark.sql.functions import row_number, monotonically_increasing_id
import pyspark.sql.functions as F


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

    @staticmethod
    def _create_argument_parser():
        """
        Create an argument parser.

        :returns
        --------
        parser - ArgumentParser
        """
        parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
            description="Python script that create an attribute list for a graph",
            epilog="Contact: stijnarend@live.nl")

        parser.version = __version__

        parser.add_argument("-d", '--data_dir', action="store",
                    dest="d", required=True, type=str,
                    help="Location of the parsed pubmed data in json format.")

        parser.add_argument("-o", '--out-dir', action="store",
                    dest="o", required=False, type=str,
                    help="The output directory to store the results - caution: requires a lot of space.")

        parser.add_argument('-v',
            '--version',
            help='Displays the version number of the script and exitst',
            action='version')

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


def create_pyspark_df(file_path: str):
    """
    Create a spark dataframe using the json files
    of the parsed pubmed articles.

    :parameters
    -----------
    file_path - str
        Path to json files of parsed articles.
    
    :returns
    --------
    df - pyspark.sql.datafram.DataFrame
        All parsed articles inside a pyspark dataframe.
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
        
    df = df.withColumn('row_id',F.monotonically_increasing_id())
    return df, sc, spark


def write_partitions_out(spark_df, partition_by: str, out_dir: str) -> None:
    """
    Write out partition of data frame to json files.

    :parameters
    -----------
    spark_df - pyspark.sql.dataframe.DataFrame
        Dataframe with the parsed article info
    partition_by - str
        Name of the column to use for partitioning
    out_dir - str
        Location of the output folder.

    Source: https://stackoverflow.com/questions/53925954/pyspark-create-multiple-json-files-from-dataframe
    """
    if not partition_by in spark_df.columns:
        print(f"Could not find column: '{partition_by}' inside supplied dataframe.")
        return None
    spark_df.write.partitionBy(partition_by).json(out_dir)
 

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
    """main"""
    cla_parser = ArgumentParser()

    data_dir = Path(cla_parser.get_argument('d'))
    out_dir = Path(cla_parser.get_argument('o'))
    # make_output_dir(out_dir)
    # json_dir = "/commons/dsls/dsph/2022/final_parsed_articles/"
    json_files = data_dir.__str__() + "/*.json"

    print("Reading data:\n")
    spark_df, sc, spark = create_pyspark_df(json_files)

    print("Writing out partitions of data frame, this wil take a while...")
    write_partitions_out(spark_df, "pmid", out_dir.__str__())

    sc.stop()


if __name__ == "__main__":
    main()
