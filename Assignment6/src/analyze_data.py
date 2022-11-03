#!/usr/bin/env python

"""Module"""

__author__ = 'Stijn Arends'
__version__ = 'v.01'

# IMPORTS
# IMPORTS
from itertools import combinations

import pandas as pd
# from fuzzywuzzy import fuzz

import pyspark
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StringType, ArrayType, DateType

import pyspark.sql.functions as F
import pyspark.sql.types as T


class AnalyzeData:

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



def main():
    """main"""
    json_dir = "/commons/dsls/dsph/2022/final_parsed_articles/"
    json_files = json_dir + "*.json"

    analyze_data = AnalyzeData(json_files)
    # print(analyze_data.schema)


if __name__ == "__main__":
    main()
