#!/usr/bin/env python3

"""
Module to parse XML files containing the PubMed literature.


To-do:
    Add try except on everything.
"""

#IMPORTS
import glob
from pathlib import Path
import re
from bs4 import BeautifulSoup
from Bio import Entrez, Medline
from numpy import str0


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
    :get_pmc
    
    To-do:
        - pass the article to the init.
        - Try except checks for empty returns
        - Add checks if the article contains the correct keys inside the init
        - try to use asynchio to read in the file
            Reading files with aynschio: https://www.twilio.com/blog/working-with-files-asynchronously-in-python-using-aiofiles-and-asyncio
            Asynchio walkthorugh: https://realpython.com/async-io-python/
        - turn the get function into asyncio functions
    """

    def __init__(self) -> None:
        pass

    def get_title(self, article) -> str:
        """
        Get the title.
        """
        return article["MedlineCitation"]['Article']['ArticleTitle']

    def get_pmid(self, article) -> str:
        """
        Get the pubmed ID.
        """
        if "PMID" in article["MedlineCitation"]:
            pmid = article["MedlineCitation"]['PMID']
        else:
            pmid = self._search_pmid_articleIdList(article['PubmedData']['ArticleIdList'] )
        return pmid

    def _search_pmid_articleIdList(self, article_id_list) -> str:
        """
        Search for the PMID inside of the ArticleIdList which is a key inside
        of the PubmedData.
        """
        pattern = re.compile("^\d+$")
        matches = [str(s) for s in article_id_list if pattern.match(s)]
        pmid = matches[0] if matches else "" # Grab the match if it found it
        return pmid

    def get_author(self, article) -> str:
        """
        Get the main author of the article.
        """
        last_name = article['MedlineCitation']['Article']['AuthorList'][0]["LastName"]
        fore_name= article['MedlineCitation']['Article']['AuthorList'][0]["ForeName"]

        return str(fore_name + " " + last_name)

    def get_co_authors(self, article) -> list[str]:
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
        return co_author_names

    def get_journal(self, article) -> str:
        """
        Get the journal of an article.
        """
        return article["MedlineCitation"]['Article']['Journal']['Title']

    def get_keywords(self, article) -> list[str]:
        """
        Get the key words of an article.
        """
        key_words = [str(word) for word in article["MedlineCitation"]['KeywordList'][0]]
        return key_words

    def get_language(self, article):
        """
        Get the langauge of an article.
        """
        return article["MedlineCitation"]['Article']['Language'][0]

    def get_doi(self, article) -> str:
        """
        Get the DOI of an article.
        """ 
        elocation_id = article["MedlineCitation"]['Article']['ELocationID']
        doi_list = elocation_id if len(elocation_id) > 0 else article['PubmedData']['ArticleIdList']
        doi = self._search_doi(doi_list)
        return doi

    def _search_doi(self, doi_list) -> str:
        """
        Search for the DOI inside of a list of potential DOIs
        using regex.

        optional pattern: '\b(10[.][0-9]{4,}(?:[.][0-9]+)*/(?:(?!["&\'<>])[[:graph:]])+)\b'
        pattern found here: https://stackoverflow.com/questions/27910/finding-a-doi-in-a-document-or-page
        """
        # Use a regex pattern to extract the DOI from the list of strings
        pattern = re.compile('(10.(\d)+\/(\S)+)')
        matches = [str(s) for s in doi_list if pattern.match(s)]
        doi = matches[0] if matches else "" # Grab the match if it found it
        return doi

    def get_pmc(self, article) -> str:
        """
        Get the PMC of the article.
        """
        article_id_list = article['PubmedData']['ArticleIdList']
        pattern = re.compile('^PMC\d+$')
        matches = [str(s) for s in article_id_list if pattern.match(s)]
        pmc = matches[0] if matches else "" # Grab the match if it found it
        return pmc


class PubmedParser:
    """
    A class for parsing pubmed articles in XML format.

    :get_publish_date
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
        # file 47:48 has elocationID
        # file 24:25 has a PMC
        for article in self.data[48:49]:
            # print(article)
            pmid = self.medline_parser.get_pmid(article)
            print(f"\nPubmed ID: {pmid}")
            title = self.medline_parser.get_title(article)
            print(f"Title: {title}")
            author = self.medline_parser.get_author(article)
            print(f"Author: {author}")
            co_authors = self.medline_parser.get_co_authors(article)
            print(f"Co-authors: {co_authors}")
            journal = self.medline_parser.get_journal(article)
            print(f"Journal: {journal}")
            key_words = self.medline_parser.get_keywords(article)
            print(f"Key words: {key_words}")
            language = self.medline_parser.get_language(article)
            print(f"Language: {language}")
            doi = self.medline_parser.get_doi(article)
            print(f"DOI: {doi}")
            pmc = self.medline_parser.get_pmc(article)
            print(f"PMC: {pmc}")


def main():
    data_dir = Path("/data/dataprocessing/NCBI/PubMed/")
    # file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0151.xml"
    # This contains one option for referecing using PMID
    file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0455.xml"

    # This contains the other option for referencing using AUthor ID
    file_2 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # Contains empty reference
    file_3 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # ---- test  1 -----

    medline_parser = MedlineParser()
    parser = PubmedParser(file, medline_parser)

    records = parser.get_parsed_xml()

    parser.parse_articles()


if __name__ == "__main__":
    main()
