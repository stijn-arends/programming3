#!/usr/bin/env python3

"""
Module to parse XML files containing the PubMed literature.
"""

#IMPORTS
import glob
from pathlib import Path
from bs4 import BeautifulSoup
from Bio import Entrez, Medline


Entrez.email = "stijnarends@live.nl"
Entrez.api_key = '9f94f8d674e1918a47cfa8afc303838b0408'


class PubmedParser:
    """
    A class for parsing pubmed articles in XML format.
    """

    def __init__(self, xml_file: Path) -> None:
        """
        Initializer
        
        :parameters
        -----------
        xml_file - Path
            XML file containing pubmed articles
        """
        self.data = self.read_pubmed_xml(xml_file)

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



def main():
    data_dir = Path("/data/dataprocessing/NCBI/PubMed/")
    # file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0151.xml"
    file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0455.xml"

    # ---- test  1 -----

    parser = PubmedParser(file)
    records = parser.get_parsed_xml()

if __name__ == "__main__":
    main()
    