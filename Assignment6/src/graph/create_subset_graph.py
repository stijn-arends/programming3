#!/usr/bin/env python

"""Module"""

__author__ = 'Stijn Arends'
__version__ = 'v.01'

# IMPORTS
import time
import pickle
import os
import sys
from typing import Any
import argparse
from pathlib import Path
import networkx as nx

import pyspark
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StructField, StringType, ArrayType, DateType

import pyspark.sql.functions as F
import pyspark.sql.types as T

SCHEMA = StructType() \
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
            description="Python script creates graphs.",
            epilog="Contact: stijnarend@live.nl")

        parser.version = __version__
        parser.add_argument("-d", '--data_dir', action="store",
                    dest="d", required=True, type=str,
                    help="Location of the parsed pubmed data in json format.")

        parser.add_argument("-n", action="store",
                           dest="n", required=False, type=int, default=2000,
                           help="Number of rows to take from the entire data frame.")

        parser.add_argument('-o', '--output', action="store",
                    dest="o", required=True,
                    help="Location of output directory.")

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
    conf = pyspark.SparkConf().setAll([('spark.executor.memory', '128g'),
                                ('spark.master', 'local[16]'),
                                ('spark.driver.memory', '128g')])
    sc = SparkContext(conf=conf)
    sc.getConf().getAll()
    spark = SparkSession(sc)
    df = spark.read.option("multiline","true").schema(SCHEMA) \
        .json(file_path)
        
    return df, sc, spark


def create_adjacency_list(x, max_size):
    """
    Helper function to create an adjacency list from
    a pyspark dataframe.

    :parameters
    -----------
    x - row
        A row of a dataframe
    max_size - int
        The maximum number of references used

    :returns
    --------
    result - Tuple[str]
        Tuple of PMIDs and None values.
    """
    id = x[0]
    vals = [None for _ in range(max_size)]
    for i in range(len(x['ref_ids'])):
        vals[i] = x['ref_ids'][i]

    return (id, *vals)


def get_data_all_articles(spark_df, adjlist_random_pd, sc):
    """
    Get the data for all the PMIDs from the pyspark dataframe this includes
    the orignal PMID that were acquired by taking the random sample of the entire dataframe
    and the PMIDs that were referenced to (available in column 'ref_ids').

    :parameters
    -----------
    spark_df - spark DataFrame
        Data frame containing info about parsed articles
    adjlist_random_pd - pd.DataFrame
        Adjacency list of the random sampled data frame
    sc - SparkContext

    :returns
    --------
    aritcles_interest_pd - pd.DataFrame
        The data of all the PMID from the random sampled df.
    """
    all_articles = []
    for i in range(len(adjlist_random_pd.columns)):
        articles = adjlist_random_pd[adjlist_random_pd.iloc[:, i].notnull()].iloc[:, i]
        all_articles.extend(articles.values)

    article_ids_broadcast = sc.broadcast(all_articles)
    articles_of_interest = spark_df.filter(F.col('pmid').isin(article_ids_broadcast.value))
    articles_interest_pd = articles_of_interest.toPandas()
    articles_interest_pd.drop_duplicates(subset=["pmid"], keep='last', inplace=True)

    return articles_interest_pd

def to_pickle(data, out_file) -> None:
    """
    Write data to a pickle file.

    :parameters
    -----------
    data - Any
        Any python object that can be pickled
    out_file - Path
        Location of output file
    """
    with open(out_file, 'wb') as file_handler:
        pickle.dump(data, file_handler)

def main():
    """main"""
    start_time = time.time()
    cla_parser = ArgumentParser()
    data_dir = cla_parser.get_argument('d')
    n_samples = cla_parser.get_argument('n')
    out_dir = Path(cla_parser.get_argument('o'))
    out_adjlist = out_dir / "adjlist_random_sample.csv"
    out_graph = out_dir / "citation_subgraph.pkl"


    # json_dir = "/commons/dsls/dsph/2022/final_parsed_articles/"
    # json_files = json_dir + "*.json"
    json_files = data_dir + "*.json"

    print("Reading data:\n")
    spark_df, sc, spark = create_pyspark_df(json_files)

    # Exract data with PMID as references
    pmid_refs = spark_df.filter(spark_df.ref_type == "pmid")

    # Take a random sample from the dataframe
    random_df = spark.createDataFrame(pmid_refs.rdd.takeSample(False, n_samples, seed=0), schema=SCHEMA)

    # Calculate the max number of references
    max_size_random = random_df.select(F.size('ref_ids')).agg({'size(ref_ids)': 'max'}).take(1)[0][0]

    # Create adjacency list from the random sampled data frame
    scheme_adjlist = StructType([StructField(f"col{i+1}", StringType(), True) for i in range(max_size_random + 1)])
    adjlist_random = random_df.rdd.map(lambda x: create_adjacency_list(x, max_size_random)).toDF(schema=scheme_adjlist)

    adjlist_random_pd = adjlist_random.toPandas()

    # Get the data from the referenced PMID as well as the source PMIDs.
    articles_interest_pd = get_data_all_articles(spark_df, adjlist_random_pd, sc)

    # Create a dictionary of node attributes
    node_attributes = articles_interest_pd.set_index('pmid').to_dict('index')

    # out_adjlist = "/commons/dsls/dsph/2022/graph_data/adjlist_random_sample.csv"

    # Write out adjacency list and read it in again to create a graph
    adjlist_random_pd.to_csv(out_adjlist, sep=' ', header=False, index=False)

    graph = nx.read_adjlist(out_adjlist, create_using=nx.DiGraph()) 

    nx.set_node_attributes(graph, node_attributes)
    # out_subset_graph = "/commons/dsls/dsph/2022/graph_data/graph_random_sample.pkl"
    to_pickle(graph, out_graph)

    elapsed_time = time.time() - start_time
    days = 0
    if elapsed_time >= 86400:
        days = int(elapsed_time / 86400)
    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print(f"\n---- Time: {days}:{elapsed} ----")


if __name__ == "__main__":
    main() # python create_subset_graph.py -d /commons/dsls/dsph/2022/final_parsed_articles/ -n 2000 -o /commons/dsls/dsph/2022/graph_data/
