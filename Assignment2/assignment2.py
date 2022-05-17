import sys
from Bio import Entrez
import multiprocessing as mp
import argparse
import os
from pathlib import Path
from typing import Any
from bs4 import BeautifulSoup

from arg_parser import ArgumentParser


__author__ = "Stijn Arends"
__version__ = "v0.1"
__data__ = "14-5-2022"


Entrez.email = "stijnarends@live.nl"
Entrez.api_key = '9f94f8d674e1918a47cfa8afc303838b0408'


class ArticleNotFound(Exception):
    """Exception raised for 404 errors.
    Attributes:
        url -- url of a website
        message -- explanation of the error
    """

    def __init__(self, pmid: int) -> None:
        self.message = f"pmid: {pmid}, is not a valid ID or it does not contain any references."
        super().__init__(self.message)


class DownloadAuthorList:
    """
    Download a number of papers that are referenced in a pubmed article. 
    """

    def __init__(self, n_articles:int, out_path: Path) -> None:
        self.n_articles = n_articles
        self.out_path = out_path

    def get_id_references(self, pmid) -> list:
        """
        Get the pubmed IDs of the references from article in pubmed.

        :returns
        --------
        references - list
            List of pubmed IDs
        """
        results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                db="pmc",
                                LinkName="pubmed_pmc_refs",
                                id=pmid))
                                
        if not results[0]["LinkSetDb"]:
            raise ArticleNotFound(pmid)

        references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]

        return references

    def download_paper(self, pmid: int) -> None:
        """
        Download an article given a pubmed id in XML format.

        :parameters
        -----------
        pmid - int
            Pubmed ID
        """
        print(f"Downloading paper: {pmid}")
        paper = Entrez.efetch(db="pmc", id=pmid, rettype="XML", retmode="text").read()
        self.get_author_names(paper)


    def get_author_names(self, paper:bytes) -> None:
        """
        Get the names of the authors of a paper.

        :parameter
        ----------
        paper - bytes
            Paper in xml format stored in a bytes object
        """
        bs4_data = BeautifulSoup(paper, 'xml')

        authors = []
        for author_info in bs4_data.find_all('name'): 
            name = " ".join(child.text for child in author_info.findChildren())
            authors.append(name)

        print(f"Authors: {authors}")


def main():
    # Get passed arguments
    cla_parser = ArgumentParser()

if __name__ == "__main__":
    main() # 30049270