#!/usr/bin/env python

"""Module"""

__author__ = 'Stijn Arends'
__version__ = 'v.01'

# IMPORTS
import pickle
from itertools import combinations
from typing import Any, Tuple
from pathlib import Path
import numpy as np
from numpy.typing import ArrayLike
import networkx as nx

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
    

    def __init__(self, file_path: str) -> None:
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
        spark_df - pyspark DataFrame
            A pyspark dataframe
        spark_context - SparkContext
        spark - SparkSession
        """
        conf = pyspark.SparkConf().setAll([('spark.executor.memory', '128g'),
                                    ('spark.master', 'local[16]'),
                                    ('spark.driver.memory', '128g')])
        spark_context = SparkContext(conf=conf)
        spark_context.getConf().getAll()
        spark = SparkSession(spark_context)
        spark_df = spark.read.option("multiline","true").schema(schema) \
            .json(file_path)
        return spark_df, spark_context, spark

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
    """
    Add doc string
    """

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

        in_degrees_sorted = in_degrees[np.argsort(in_degrees)[::-1]]
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

    def author_similar_authors_refs(self) -> Tuple[list, dict, float]:
        """
        Find out if authors mainly reference papers with other authors with whom
        they've co-authored papers (including themselves).

        :returns
        --------
        overlap_percentages - list
            List of overlapping percentages
        highest_perc_overlap - dict
            Author with the highest % overlap
        avg_overlap - float
            The average % overlap
        """
        author_co_authors_flat = self._get_author_co_authors_flat()
        authors_with_ref = self.get_authors_with_ref()

        author_co_auth_refs = {author:{'total': 0, 'yes': 0} for author in authors_with_ref.values()}

        for node, author in authors_with_ref.items():
            out_degree = self.graph.out_degree(node)
            co_auths = author_co_authors_flat[author]

            if out_degree != 0:
                for neighbour in self.graph.neighbors(node):
                    co_auths_refs = self.graph.nodes[neighbour]['co_authors']
                    co_auths_refs.insert(0, self.graph.nodes[neighbour]['author'])
                    author_co_auth_refs[author]['total'] += 1
                    if list(set(co_auths_refs) & set(co_auths)):
                        author_co_auth_refs[author]['yes'] += 1

                        
        overlap_percentages, highest_perc, avg_overlap = self._get_overlap_percentages(author_co_auth_refs)
        return overlap_percentages, highest_perc, avg_overlap

    def _get_author_co_authors_flat(self) -> dict[list]:
        """
        Get all the co authors of all the articles for each other.
        Store these in one lists.

        :returns
        --------
        author_co_authors_flat - dict
            authors and there co authors for each of their publications
        """
        author_co_authors_flat = {}

        for pmid, author in self.all_authors.items():
            co_authors = self.graph.nodes[pmid]['co_authors']
            if author not in author_co_authors_flat:
                author_co_authors_flat[author] = co_authors
                author_co_authors_flat[author].insert(0, author)
            else:
                author_co_authors_flat[author].extend(co_authors)
                # author_co_authors_flat[author].insert(0, author)

        return author_co_authors_flat

    def get_authors_with_ref(self) -> dict[str]:
        """
        Get the authors and articles that have references.

        :returns
        --------
        authors_with_ref - dict
            PMID and author name of articles that contain references
        """
        authors_with_ref = {}
        for article, author in self.all_authors.items():
            out_degr = self.graph.out_degree(article)
            if out_degr != 0:
                authors_with_ref[article] = author

        return authors_with_ref

    @staticmethod
    def _get_overlap_percentages(author_co_auth_refs) -> Tuple[list, dict, float]:
        """
        Get the percentage of overlapping authors/co-authors from references.

        :parameters
        -----------
        author_co_auth_refs - dict
            dict of authors with stats on the number of articles they have
            and the number of times those articles had author/co-authors
            that they often work with

        :returns
        --------
        overlap_percentages - list
            List of overlapping percentages
        highest_perc_overlap - dict
            Author with the highest % overlap
        avg_overlap - float
            The average % overlap
        """
        highest_perc_overlap = {'author': None, 'overlap': 0, 'total': 0}

        overlap_percentages = []

        for author, stats in author_co_auth_refs.items():
            overlap = (stats['yes'] / stats['total']) * 100
            overlap_percentages.append(overlap)
            if overlap > highest_perc_overlap['overlap']:
                highest_perc_overlap['overlap'] = overlap
                highest_perc_overlap['author'] = author
                highest_perc_overlap['total'] = stats['total']

        avg_overlap = np.mean(overlap_percentages)
        print(f"Mean overlap: {avg_overlap:.3}%")
        return overlap_percentages, highest_perc_overlap, avg_overlap

    def distinguish_citations(self, threshold:float = 1) -> Tuple[ArrayLike, ArrayLike]:
        """
        Distinguish what are highly and lowly cited papers.

        :parameters
        -----------
        threshold - float
            The threshold to decide what is high and what is low,
            default = 1

        :returns
        --------
        highly_cited - ArrayLike
            highly cited articles
        lowly_cited - ArrayLike
            lowly cited articles
        """
        in_degrees = np.array(list(dict(self.graph.in_degree).values()))
        articles = np.array(list(dict(self.graph.in_degree).keys()))

        highly_cited = articles[np.where(in_degrees >= threshold)[0]]
        lowly_cited = articles[np.where(in_degrees < threshold)[0]]
        return highly_cited, lowly_cited


    def get_time_span(self, articles: ArrayLike) -> list:
        """
        Get the time span of reference for a list of articles.

        :parameters
        -----------
        articles - ArrayLike
            Array of articles (nodes)

        :returns
        --------
        time_span - list
            Time span of being referenced for each article
        """
        time_span = []

        for article in articles:
            all_years = []
            in_edges = self.graph.in_edges(article, data=True)
            if in_edges:
                for ref, _, _ in in_edges:
                    year = self.graph.nodes[ref]['publish_date'].year
                    all_years.append(year)
                
                span = np.max(all_years) - np.min(all_years)
                time_span.append(span)
        return time_span

    @staticmethod
    def find_intersect_out_neighbours(graph: nx.DiGraph, node: str, attribute: str,
            total: int, intersect: int) -> Tuple[int, int]:
        """
        Find if there are similarities for the specified attributes between the source node and
        the out going neighbours (i.e. articles that have been cited by the source paper).

        :parameters
        -----------
        graph - nx.DiGraph
            A citation network
        node - str
            Name of node
        attribute - str
            Name of node attribute
        total - int
            Total number of articles parsed
        intersect - int
            Total number of intersect found
        """
        attribute_source = graph.nodes[node][attribute]
        for out_neighbor in graph.neighbors(node):
            attribute_ref = graph.nodes[out_neighbor][attribute]
            if not attribute_ref:
                continue
            total += 1
            # print(key_words_out_neigh)
            if list(set(attribute_source) & set(attribute_ref)):
                intersect += 1
        return total, intersect

    @staticmethod
    def find_intersect_in_neighbours(graph: nx.DiGraph, node: str, attribute: str,
            total: int, intersect: int) -> Tuple[int, int]:
        """
        Find if there are similarities for the specified attributes between the source node and
        the in going neighbours (i.e. articles that have cited the source paper).

        :parameters
        -----------
        graph - nx.DiGraph
            A citation network
        node - str
            Name of node
        attribute - str
            Name of node attribute
        total - int
            Total number of articles parsed
        intersect - int
            Total number of intersect found
        """
        attribute_source = graph.nodes[node][attribute]
        for ref, _, _ in graph.in_edges(node, data=True):
            attribute_ref = graph.nodes[ref][attribute]
            if not attribute_ref:
                continue

            total += 1
            if list(set(attribute_source) & set(attribute_ref)):
                intersect += 1
        return total, intersect

    def find_corr_citation_key_words(self, nodes: ArrayLike) -> list[float]:
        """
        Find if there is a correlation between the citations and the number of key words that paper share.
        I.e. papers which share the same subject cite each other more often.

        :parameters
        -----------
        nodes - ArrayLike
            An array of nodes

        :returns
        --------
        correlation_citaion_key_words - list
            The correlation between citations and number of key words
            for each article (node)
        """
        if isinstance(nodes, str):
            nodes = [nodes]

        correlation_citation_key_words = []
        for node in nodes:
            out_degree = self.graph.out_degree(node)
            in_degree = self.graph.in_degree(node)

            # If there are no references skip
            if out_degree == 0 and in_degree == 0:
                continue

            # If there are no key words skip
            key_words = self.graph.nodes[node]['key_words']
            if not key_words:
                continue
            total = 0
            intersect = 0

            # Check if the outgoing neighbours share similar key words
            if out_degree != 0:
                # print('neighbors:')
                total, intersect = self.find_intersect_out_neighbours(self.graph, node, 'key_words', total, intersect)

            # Check if ingoing neighbours share similar key words
            if in_degree != 0:
                total, intersect = self.find_intersect_in_neighbours(self.graph, node, 'key_words', total, intersect)
            if total != 0:
                correlation = intersect / total
                correlation_citation_key_words.append(correlation)
        return correlation_citation_key_words

    def get_author_citations(self) -> dict:
        """
        Get for each article the name of the author and the number
        of citations it has.

        :returns
        --------
        author_citations - dict
            Author name and number of citations it got for an article.
        """
        author_citations = {"author": [], 'citations': []}

        for article, author in self.all_authors.items():
            in_degree = self.graph.in_degree(article)
            if author != '':
                author_citations["author"].append(author)
                author_citations['citations'].append(in_degree)

        author_citations_df = pd.DataFrame(author_citations)
        return author_citations_df

    def most_cited_author(self) -> Tuple[str, int]:
        """
        Get the name of the most cited author as well as the number of
        citations that they got.

        :returns
        --------
        author - str
            Name of the most cited author
        max_n_citations - int
            Number of citations that author got
        """
        author_citations_df = self.get_author_citations()
        grouped_author_citation = author_citations_df.groupby('author').agg({'citations':'sum'})
        max_n_citations = grouped_author_citation.citations.max() 
        most_cited_author = grouped_author_citation[grouped_author_citation['citations'] == max_n_citations]
        return most_cited_author.index[0], max_n_citations

    def calcualte_h_index(self) -> Tuple[list, int]:
        """
        Calculate the h-index for each author and return the name of
        the author(s) that got the highest h-index.

        The h-index is a measure of the number of publications published (productivity),
        as well as how often they are cited.

        h-index = the number of publications with a citation number greater than or equal to h.

        For example:
            * 15 publications cited 15 times or more, is a h-index of 15.

        :returns
        --------
        authors - list
            List of author(s) with the highest h-index
        max_h_index - int
            The highest h-index
        """
        author_citations_df = self.get_author_citations()
        author_citations_df['h-index'] = author_citations_df.groupby('author')['citations']\
            .transform( lambda x: (x >= x.count()).sum())

        max_h_index = author_citations_df['h-index'].max()
        max_h_index_df = author_citations_df[author_citations_df['h-index'] == max_h_index]

        authors = max_h_index_df.author.unique().tolist()
        return authors, max_h_index

    def find_relation_language_references(self) -> Tuple[list, float]:
        """
        Find how often papers reference other papers that were written in
        the same language.

        :returns
        --------
        overlap_languages - list
            The amount of references that were written using the same
            language depicted in %.
        avg_overlap_languages - float
            The average percentage
        """
        overlap_languages = []
        
        for node in self.graph.nodes():
            out_degr = self.graph.out_degree(node)
            language_source = self.graph.nodes[node]['language']
            total, intersect = 0, 0
            if out_degr != 0:
                for neighbor in self.graph.neighbors(node):
                    total += 1
                    language_neighbor = self.graph.nodes[neighbor]['language']
                    if language_source == language_neighbor:
                        intersect += 1
        
                overlap = (intersect / total) * 100
                overlap_languages.append(overlap)
        return overlap_languages, np.mean(overlap_languages)

def make_data_dir(path: Path) -> None:
    """
    Create a directory (if it does not exsit yet) to store the
    data.

    :Excepts
    --------
    FileExistsError
        The directory already exists
    """
    try:
        path.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(f"[{make_data_dir.__name__}] {path} already exists.")


def main():
    """main"""
    json_dir = "/commons/dsls/dsph/2022/final_parsed_articles/"
    json_files = json_dir + "*.json"
    subset_graph = "/commons/dsls/dsph/2022/graph_data/citation_subgraph_3.pkl"

    out_path = Path(__file__).parent.parent.absolute()
    result_files = out_path / 'result_files'
    make_data_dir(result_files)
    print(f"Output path: {out_path}") 

    analyze_data = AnalyzeData(json_files)

    # Q1
    avg_co_authors = analyze_data.avg_co_authors()
    print(f"Average co authors: {avg_co_authors}")

    graph = read_pickl(subset_graph)
    analyze_graph = AnalyzeGraph(graph)
    most_cited_paper, max_citation_paper, in_degrees_sorted = analyze_graph.most_cited_paper(additional_info=True)

    threshold = np.percentile(in_degrees_sorted, 99.5)
    print(f"Threshold: {threshold}")

    print(f"Most cited paper: {most_cited_paper}")
    # Q2
    relation_auth_co_auths, avg_relation_co_auths = analyze_graph.author_similar_co_authors()
    # Q3
    overlap_percentages, highest_perc, avg_overlap = analyze_graph.author_similar_authors_refs()
    # Q4
    highly_cited, lowly_cited = analyze_graph.distinguish_citations(threshold=threshold)
    # Save this
    highly_cited_time_span = analyze_graph.get_time_span(highly_cited)
    np.save(result_files / 'highly_cited_time_span.npy', highly_cited_time_span)

    avg_time_span_high = np.mean(highly_cited_time_span)
    print(f"Average time span highly cited papers: {avg_time_span_high:.4f} years")

    # Save this
    lowly_cited_time_span = analyze_graph.get_time_span(lowly_cited)
    np.save(result_files / 'lowly_cited_time_span.npy', lowly_cited_time_span)
    avg_time_span_low = np.mean(lowly_cited_time_span)
    print(f"Average time span lowly cited papers: {avg_time_span_low:.4f} years")

    # Save this:
    time_span_highest = [graph.nodes[ref]['publish_date'].year \
        for ref, _, _ in graph.in_edges(most_cited_paper, data=True)]
    np.save(result_files / 'time_span_most_cited_paper.npy', time_span_highest)
    

    # Q5
    correlation_citation_key_words = analyze_graph.find_corr_citation_key_words(graph.nodes())
    avg_corr_citation_key_words = np.mean(correlation_citation_key_words)
    print(f"Mean corr: {avg_corr_citation_key_words:.3f}")
    print(f"Number of acceptable articles: {len(correlation_citation_key_words)}")

    # Q6
    corr_citation_key_words_highly_cited = analyze_graph.find_corr_citation_key_words(highly_cited)
    avg_corr_citation_key_words_high = np.mean(corr_citation_key_words_highly_cited)
    print(f"Mean corr: {np.mean(corr_citation_key_words_highly_cited):.3f}")

    # Q7: h index
    authors, max_h_index = analyze_graph.calcualte_h_index()
    print(f"Author(s) with highest h-index: {authors} - {max_h_index}")

    # Q8: most cited author
    most_cited_author, max_citation_author = analyze_graph.most_cited_author()
    print(f'Most cited author: {most_cited_author}')

    # Q9: Do papers mostly reference other papers from the same language?
    overlap_languages, avg_overlap_languages = analyze_graph.find_relation_language_references()


    questions = ["How large a group of co-authors does the average publication have?",
                "Do authors mostly publish using always the same group of authors?",
                "Do authors mainly reference papers with other authors with whom they've " \
                    "co-authored papers (including themselves)?",
                "What is the distribution in time for citations of papers in general, and "\
                    "for papers with the highest number of citations? Do they differ?",
                "Is there a correlation between citations and the number of keywords that "\
                    "papers share? I.e. papers which share the same subject cite each other more often.",
                "For the most-cited papers (define your own cutoff), is the correlation in "\
                    "shared keywords between them and the papers that cite them different from (5) ?",
                "Do papers mostly reference other papers that were written in the same language?",
                "What is the most cited paper?", "Who is the most cited author?",
                "Which author(s) has/have the highest h-inex?"]

    most_cited_paper_info = {'PMID': most_cited_paper, 'author': graph.nodes[most_cited_paper]['author'],
                'citations': max_citation_paper, 'title': graph.nodes[most_cited_paper]['title']}

    most_cited_author_info = {'author': most_cited_author, 'citations': max_citation_author}

    answers = [avg_co_authors, f"{np.round(avg_relation_co_auths, 3)}%", f"{np.round(avg_overlap, 3)}%",
            {'general': f"{np.round(avg_time_span_low, 4)} years", 'high': f"{np.round(avg_time_span_high, 2)} years"},
            np.round(avg_corr_citation_key_words, 3), np.round(avg_corr_citation_key_words_high, 3),
            f"{avg_overlap_languages:.3}%",
            most_cited_paper_info, most_cited_author_info, {'authors': authors, 'h-index': max_h_index}]

    questions_answers = pd.DataFrame({'Question': questions, 'Asnwer': answers})
    questions_answers.to_csv(out_path / 'questions_answers.csv', sep=',', header=True)



if __name__ == "__main__":
    main()
