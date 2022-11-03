#!/usr/bin/env python

"""Module"""

__author__ = 'Stijn Arends'
__version__ = 'v.01'

# IMPORTS
import pickle
from itertools import combinations
from typing import Any, Tuple
import numpy as np
from numpy.typing import ArrayLike
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
# from fuzzywuzzy import fuzz

import pyspark
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StringType, ArrayType, DateType

import pyspark.sql.functions as F
import pyspark.sql.types as T


def read_pickl(file) -> Any:
    """
    Read in a pickle file.
    """
    with open(file, 'rb') as fh:
        data = pickle.load(fh)

    return data


class AnalyzeData:
    """
    add docstring
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
    

    def __init__(self, file_path) -> None:
        self.file_path = file_path
        self.spark_df, self.sc, self.spark = self.read_spark_df_json(file_path, self.schema)

    @staticmethod
    def read_spark_df_json(file_path: str, schema: StructType):
        """
        Read in a pyspark data frame from multiple json file.

        :parameters
        -----------
        file_path - str
            Directory containg json files used to create dataframe
        schema - StructType
            Schema to infer for the data

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
        df = spark.read.option("multiline","true").schema(schema) \
            .json(file_path)
        return df, sc, spark

    def avg_co_authors(self) -> float:
        """
        Calculate the average amount of co-authors each article has.

        :returns
        --------
        avg_co_authors - float
            Average amount of co-authors per publication
        """
        count_df = self.spark_df.select('co_authors',F.size('co_authors').alias('n_co_authors'))
        avg_co_authors = count_df.select(F.mean('n_co_authors')).take(1)[0][0]
        return avg_co_authors


class AnalyzeGraph:

    def __init__(self, graph) -> None:
        self.graph = graph
        self.all_authors = nx.get_node_attributes(graph, 'author')

    def most_cited_paper(self, additional_info: bool = False) -> Tuple[str, int, ArrayLike]:
        """
        Get the most cited paper. If additional_info is set to True
        (default False) then the function will also return the number
        of citations and a list of articles sorted on number of citations
        from highest to lowest.

        :parameters
        -----------
        additional_info - bool
            Wheter or not to return the additional info.

        :returns
        --------
        most_cited_paper - str
            PMID of the most cited paper
        max_citation - int
            Number of maximum citations a paper got
        in_degrees_sorted - np.array
            List of articles sorted based on number of citations from highest
            to lowest
        """
        in_degrees = np.array(list(dict(self.graph.in_degree).values()))
        index_max = np.argmax(in_degrees)
        most_cited_paper, max_citation = list(self.graph.in_degree)[index_max]
        if not additional_info:
            return most_cited_paper

        in_degrees_sorted = np.argsort(in_degrees)[::-1]
        return most_cited_paper, max_citation, in_degrees_sorted

    def author_similar_co_authors(self) -> Tuple[dict, float]:
        """
        Find out if authors mainly puplish using the same group of
        co-authors by calculating the amount of times there were similarities
        in the co-author list for different publications.

        :returns
        --------
        relation_auth_co_auths - dict
            the frequency of authors publishing with similar co authors
        avg_relation_co_auths - float
            average frequency
        """
        author_co_authors = self._get_author_co_authors()
        relation_auth_co_auths = {}

        for author, co_authors_list in author_co_authors.items():
            if len(co_authors_list) > 1:
                comb = list(combinations(co_authors_list, 2))
                total = len(comb)
                intersect = 0
                for list_one, list_two in comb:
                    # use this to check if there is an intersect: 
                    if list(set(list_one) & set(list_two)):
                        intersect += 1
                percentage = (intersect / total) * 100
                relation_auth_co_auths[author] = percentage
        avg_relation_co_auths = np.array(list(relation_auth_co_auths.values())).mean()
        print(f"Average percentage of authors publishing with the same or similar group of authors is: {avg_relation_co_auths:.3f}%")

        return relation_auth_co_auths, avg_relation_co_auths


    def _get_author_co_authors(self) -> dict[list]:
        """
        Get all the co authors of all the articles for each other.
        Store these in individual lists.

        :returns
        --------
        author_co_authors - dict
            authors and there co authors for each of their publications
        """
        author_co_authors = {}

        for pmid, author in self.all_authors.items():
            co_authors = self.graph.nodes[pmid]['co_authors']
            if author not in author_co_authors:
                author_co_authors[author] = [co_authors]
            else:
                author_co_authors[author].append(co_authors)

        return author_co_authors


def main():
    """main"""
    json_dir = "/commons/dsls/dsph/2022/final_parsed_articles/"
    json_files = json_dir + "*.json"
    subset_graph = "/commons/dsls/dsph/2022/graph_data/citation_subgraph_3.pkl"


    # analyze_data = AnalyzeData(json_files)

    graph = read_pickl(subset_graph)
    analyze_graph = AnalyzeGraph(graph)
    print(analyze_graph.most_cited_paper())
    analyze_graph.author_similar_co_authors()

if __name__ == "__main__":
    main()
