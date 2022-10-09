#!/usr/bin/env python3

"""
Module to parse XML files containing the PubMed literature.


To-do:
    Add try except on everything.
"""

#IMPORTS
import glob
from pathlib import Path
from bs4 import BeautifulSoup
from Bio import Entrez, Medline


Entrez.email = "stijnarends@live.nl"
Entrez.api_key = '9f94f8d674e1918a47cfa8afc303838b0408'


class MedlineParser:
    """
    Parse MedlineCitation of pubmed article.

    :get_title
    :get_pmid
    :get_author
    :get_co_author
    :get_language
    :get_article_date
    :get_keywords
    :
    """

    def __init__(self) -> None:
        pass

    def get_title(self, article):
        """
        Get the title.
        """
        return article["MedlineCitation"]['Article']['ArticleTitle']

    def get_pmid(self, article):
        """
        Get the pubmed ID.
        """
        return article["MedlineCitation"]['PMID']

    def get_author(self, article):
        """
        Get the main author of the article.
        """
        last_name = article['MedlineCitation']['Article']['AuthorList'][0]["LastName"]
        fore_name= article['MedlineCitation']['Article']['AuthorList'][0]["ForeName"]

        return str(fore_name + " " + last_name)

    def get_co_authors(self, article):
        """
        Get the co authors
        """
        authors = article['MedlineCitation']['Article']['AuthorList'][1:]
        co_author_names = []
        if authors:
            for author in authors:
                last_name = author["LastName"]
                fore_name= author["ForeName"]
                co_author_names.append(str(fore_name + " " + last_name))

        # print(article['MedlineCitation']['Article']['AuthorList'])
        # print(co_author_names)
        return co_author_names


class PubmedParser:
    """
    A class for parsing pubmed articles in XML format.
    """

    def __init__(self, xml_file: Path, medline_parser: MedlineParser) -> None:
        """
        Initializer
        
        :parameters
        -----------
        xml_file - Path
            XML file containing pubmed articles
        """
        self.data = self.read_pubmed_xml(xml_file)
        self.medline_parser = medline_parser

    @staticmethod
    def read_pubmed_xml(file: Path) -> list:
        """
        Read a pubmed file in xml format that contains articles.

        :parameters
        -----------
        file - Path
            XML file containing pubmed articles
        
        :returns
        --------
        records - list
            List containing information about articles inside of dictionaries.
        """
        with open(file, "br") as handle:
            records = Entrez.read(handle)

        return records['PubmedArticle']

    def get_parsed_xml(self) -> list:
        """
        Get the parsed XML file.

        :returns
        --------
        self.data - list
            Parsed XML data containing info about the articles
        """
        return self.data


    def parse_articles(self):
        """
        Create a dictionary of lists and save the results in here.
        """
        for article in self.data[:1]:
            print(article)
            pmid = self.medline_parser.get_pmid(article)
            print(f"\nPubmed ID: {pmid}")
            title = self.medline_parser.get_title(article)
            print(f"Title: {title}")
            author = self.medline_parser.get_author(article)
            print(f"Author: {author}")
            co_authors = self.medline_parser.get_co_authors(article)
            print(f"Co-authors: {co_authors}")




def main():
    data_dir = Path("/data/dataprocessing/NCBI/PubMed/")
    # file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0151.xml"
    file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0455.xml"

    # ---- test  1 -----

    medline_parser = MedlineParser()
    parser = PubmedParser(file, medline_parser)

    records = parser.get_parsed_xml()

    parser.parse_articles()


if __name__ == "__main__":
    main()
