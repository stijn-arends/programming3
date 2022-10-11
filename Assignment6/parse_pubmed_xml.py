#!/usr/bin/env python3

"""
Module to parse XML files containing the PubMed literature.


To-do:
    Add try except on everything.
"""

#IMPORTS
# import glob
from pathlib import Path
import re
import datetime
from typing import Any
# from bs4 import BeautifulSoup
from Bio import Entrez


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
            Reading files with aynschio:
            https://www.twilio.com/blog/working-with-files-asynchronously-in-python-using-aiofiles-and-asyncio
            Asynchio walkthorugh: https://realpython.com/async-io-python/
        - turn the get function into asyncio functions
        - Move this class to a seperate module
    """

    def __init__(self, article: dict) -> None:
        """
        Initilizer of the MedlineParser class.

        Check if the article contains the correct keys.
        """
        self.article = article

        assert "MedlineCitation" in article, "Key: MedLineCitation not found inside the article."
        assert "PubmedData" in article, "Key: PubmedData not found inside the article."

    def get_title(self) -> str:
        """
        Get the title.

        :returns
        --------
        pmid - str
            The title of the article, returns empty str when no title could be found
        """
        try:
            title = self.article["MedlineCitation"]['Article']['ArticleTitle']
        except IndexError:
            title = ""
        return title

    def get_pmid(self) -> str:
        """
        Get the pubmed ID. First search the 'PMID' key inside the 'MedlineCitation'
        section. If it can't be found there try to extract the PMID from the
        'ArticleIdList' inside the PubmedData using regex.

        :returns
        --------
        pmid - str
            The PMID of the article, returns empty str when no PMID could be found

        :see also
        --------
        MedlineParser._search_pmid_article_list()
        """
        try:
            if "PMID" in self.article["MedlineCitation"]:
                pmid = self.article["MedlineCitation"]['PMID']
            else:
                pmid = self._search_pmid_article_list(self.article['PubmedData']['ArticleIdList'])
        except IndexError:
            pmid = ""
        return pmid

    @staticmethod
    def _search_pmid_article_list(article_id_list) -> str:
        """
        Search for the PMID inside of the ArticleIdList which is a key inside
        of the PubmedData.
        """
        pattern = re.compile(r"^\d+$")
        matches = [str(s) for s in article_id_list if pattern.match(s)]
        pmid = matches[0] if matches else "" # Grab the match if it found it
        return pmid

    def get_author(self) -> str:
        """
        Get the main author of the article.

        :returns
        --------
        full_name - str
            Full name of the author, returns empty str when author name could not be found
        """
        try:
            last_name = self.article['MedlineCitation']['Article']['AuthorList'][0]["LastName"]
            initials = self.article['MedlineCitation']['Article']['AuthorList'][0]["Initials"]

            # add . between and after the initals
            if '.' not in initials:
                initials = ".".join(initials) + '.'
            full_name = str(last_name + ", " + initials)
        except IndexError:
            full_name = ""

        return full_name

    def get_co_authors(self) -> list[str]:
        """
        Get the co authors of the article.

        :returns
        --------
        co_authors - list[str]
            A list of the co authors, returns empty list when no co authors could be found
        """
        try:
            authors = self.article['MedlineCitation']['Article']['AuthorList'][1:]
            co_author_names = []
            if authors:
                for author in authors:
                    last_name = author["LastName"]
                    initials = author["Initials"]
                    # add . between and after the initals
                    if '.' not in initials:
                        initials = ".".join(initials) + '.'
                    co_author_names.append(str(last_name + ", " + initials))
        except IndexError:
            co_author_names = []
        return co_author_names

    def get_journal(self) -> str:
        """
        Get the journal of an article.

        :returns
        --------
        journal - str
            The journal of the article, returns empty str when journal could not be found
        """
        try:
            journal = self.article["MedlineCitation"]['Article']['Journal']['Title']
        except IndexError:
            journal = ""
        return journal

    def get_keywords(self) -> list[str]:
        """
        Get the key words of an article.

        :returns
        --------
        key_words - list[str]
            A list of the key words, returns empty list when no key words could be found
        """
        try:
            key_words = [str(word) for word in self.article["MedlineCitation"]['KeywordList'][0]]
        except IndexError:
            key_words = []
        return key_words

    def get_language(self):
        """
        Get the langauge of an article.

        :returns
        --------
        language - str
            The language that was used to write the article,
            returns empty str when language could not be found
        """
        try:
            language = self.article["MedlineCitation"]['Article']['Language'][0]
        except IndexError:
            language = ""
        return language

    def get_doi(self) -> str:
        """
        Get the DOI of an article.

        First, search inside the 'ELocationID' key check if it is not empty and if
        that is the case store the value. If it is empty look in the ArticleIdList key and
        store that value. Next, using regex extract the DOI.

        :returns
        --------
        doi - str
            the Digital Object Identifier(DOI) of the article.

        :see also
        --------
        MedlineParser._search_doi()
        """
        try:
            elocation_id = self.article["MedlineCitation"]['Article']['ELocationID']
            article_id_list = self.article['PubmedData']['ArticleIdList']
            doi_list = elocation_id if len(elocation_id) > 0 else article_id_list
            doi = self._search_doi(doi_list)
        except IndexError:
            doi = ""
        return doi

    @staticmethod
    def _search_doi(doi_list) -> str:
        """
        Search for the DOI inside of a list of potential DOIs
        using regex.

        optional pattern: '\b(10[.][0-9]{4,}(?:[.][0-9]+)*/(?:(?!["&\'<>])[[:graph:]])+)\b'
        pattern found here:
        https://stackoverflow.com/questions/27910/finding-a-doi-in-a-document-or-page

        :parameters
        -----------
        doi_list - list
            List of potential DOI strings

        :returns
        --------
        doi - str
            the Digital Object Identifier(DOI) of the article.
        """
        # Use a regex pattern to extract the DOI from the list of strings
        pattern = re.compile(r'(10.(\d)+\/(\S)+)')
        matches = [str(s) for s in doi_list if pattern.match(s)]
        doi = matches[0] if matches else "" # Grab the match if it found it
        return doi

    def get_pmc(self) -> str:
        """
        Get the Pubmed Central (PMC) identifier of the article.

        Look trough the 'ArticleIdList' and try to find the PMC identifier
        using regex. If nothing is found an empty string is returned.

        :returns
        --------
        pmc - str
            Either the PMC identifier of an article or an empty string if
            it does not have one.
        """
        try:
            article_id_list = self.article['PubmedData']['ArticleIdList']
            pattern = re.compile(r'^PMC\d+$')
            matches = [str(s) for s in article_id_list if pattern.match(s)]
            pmc = matches[0] if matches else "" # Grab the match if it found it
        except IndexError:
            pmc = ""
        return pmc


class PubmedParser:
    """
    A class for parsing pubmed articles in XML format.

    :get_publish_date
    :get_references

    To-do:
        - Think about using the __iter__ method instead of the
        parse_articles method.
        Inspiration mcoding video:
        https://www.youtube.com/watch?v=tmeKsb2Fras
        Look at prog2 exercises for more inspiration.
        - Decide which reference type you want
            - all authors
            - only main author
        Only the main author should be sufficient for coupling authors
        to the PMID after everything has been parsed. 
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

    def parse_articles(self):
        """
        Create a dictionary of lists and save the results in here.

        for file or file_1:
            * article 24:25 has a PMC
            * article 47:48 has elocationID

        For the 'file_2':
            * article 0:1 has no ReferenceList
            * article 12:13 has a citation, articleIdList reference
            * article 5403:5404 has a citation reference
        """
        # file 47:48 has elocationID
        # file 24:25 has a PMC
        for article in self.data[12:13]:
            # print(article)
            medline_parser = MedlineParser(article)
            pmid = medline_parser.get_pmid()
            print(f"\nPubmed ID: {pmid}")
            title = medline_parser.get_title()
            print(f"Title: {title}")
            author = medline_parser.get_author()
            print(f"Author: {author}")
            co_authors = medline_parser.get_co_authors()
            print(f"Co-authors: {co_authors}")
            journal = medline_parser.get_journal()
            print(f"Journal: {journal}")
            key_words = medline_parser.get_keywords()
            print(f"Key words: {key_words}")
            language = medline_parser.get_language()
            print(f"Language: {language}")
            doi = medline_parser.get_doi()
            print(f"DOI: {doi}")
            pmc = medline_parser.get_pmc()
            print(f"PMC: {pmc}")
            pubblish_date = self.get_publish_date(article)
            print(f"Publish date: {pubblish_date}")

            ref_id_list = self.get_references(article)
            print(f"Reference ID list: {ref_id_list}")

    @staticmethod
    def get_publish_date(article: dict) -> datetime.date:
        """
        Get the date of publishing the article in year-month-day format.

        :parameters
        -----------
        article - dict
            The article inside a dictionary.

        :returns
        --------
        date - datetime.date
            The date of publishing, returns empty str when date can't be found.
        """
        try:
            year = int(article['PubmedData']['History'][0]['Year'])
            month = int(article['PubmedData']['History'][0]['Month'])
            day = int(article['PubmedData']['History'][0]['Day'])
            date = datetime.date(year, month, day)
        except IndexError:
            date = ""
        return date

    def get_references(self, article: dict) -> list[Any]:
        """
        Extract either the main author or the PMID of all the references.
        First look if the article contains references. Next, check which reference 
        style was used:

        1.	A dictionary with 2 keys, “citation” and “ArticleIdList”. The citation has as value:
            journal + journal edition, and ArticleIdList has as value: a string representation
            of the PMID.

            Example:
                {'Citation': 'J Biol Chem. 1958 Nov;233(5):1135-9', 'ArticleIdList':
                [StringElement('13598746', attributes={'IdType': 'pubmed'})]}

        2.	A dictionary with 1 key: “citation”, which has a value containing:
            authors + year + title + journal edition

            Example:
                {'Citation': 'Adam, E.J.H. and Adam, S.A. 1994. Identification of
                cytosolic factors required for nuclear location sequence-mediated
                binding to the nuclear envelope. J. Cell Biol. 125:547-555.'}

        :parameters
        -----------
        article - dict
            The article inside a dictionary.

        :returns
        --------
        ref_id_list - list[str]
            List of either the PMID of the references or the authors.

        :see also
        ---------
        PubmedParser._get_pmid_reference()
        PubmedParser._get_authors_reference()
        PubmedParser._get_main_author_references()
        """
        # print(article['PubmedData']['ReferenceList'])
        try:
            reference_list = article['PubmedData']['ReferenceList'][0]
            #  print(f"references: {len(reference_list)}")

            references = reference_list['Reference']
            if len(references[0].keys()) > 1:
                # print(f"\nTest: {references[0]['ArticleIdList'][0].attributes}")
                # print(len(references))
                ref_id_list = self._get_pmid_reference(references)

                # print(f"PMID list: {ref_id_list}")
                # print(len(ref_id_list))
            else:
                # ref_id_list = self._get_authors_reference(references)
                # print(f"Authors list references: {ref_id_list}\n")
                ref_id_list = self._get_main_author_references(references)
                # print(f"Main authors list references: {ref_id_list}\n")
        except IndexError:
            ref_id_list = []

        return ref_id_list

    @staticmethod
    def _get_pmid_reference(references: list) -> list[int]:
        """
        Get the PMID from the ArticleIdList for each reference.

        :parameters
        -----------
        reference - list
            List of the references.

        :returns
        --------
        pmid_list - list[str]
            List of the PMID from the references.
        """
        pmid_list = []
        for reference in references:
            # print(reference)
            # print(reference["ArticleIdList"][0])
            article_id_list = reference["ArticleIdList"][0]
            pmid = int(article_id_list) if article_id_list.attributes["IdType"] == "pubmed" else ""
            # print(f"PMID: {pmid}")
            pmid_list.append(pmid)
        
        return pmid_list

    @staticmethod
    def _get_authors_reference(references: list) -> list[str]:
        """
        Get the names of the authors for each reference using regex. Extract all the information
        before the first number is found (i.e. year of publishing). Next, split that text on '. '
        and if that did not work split it on ', ' to seperate the different authors. If this is the
        case then combine the sir name and initials and put them into a new list. 
        
        :parameters
        -----------
        reference - list
            List of the references.

        :returns
        --------
        authors_references - list[str]
            List of the authors from the references.

        Format of reference: 
            * authors + year + title + journal edition / book

        Example:
            'Adam, E.J.H. and Adam, S.A. 1994.
            Identification of cytosolic factors required for nuclear location sequence-mediated
            binding to the nuclear envelope. J. Cell Biol. 125:547-555.'
        """
        pattern = r'^(.+?)(?=\d).*'
        authors_references = []
        try: 
            for reference in references:
                check = re.search(pattern, reference["Citation"])
                names = check.group(1)
                
                names = names.replace('and ', '').rstrip()
                names = names.split('. ')
                correct_names = []
                if len(names) == 1:
                    names = names[0].split(', ')
                    for i in range(1, len(names), 2):
                        correct_names.append(names[i -1] + ", " + names[i])

                    authors_references.append(correct_names)
                else:
                    authors_references.append(names)
        except Exception:
            pass

        # print(f"Authors of the references: {authors_references}")
        return authors_references

    @staticmethod
    def _get_main_author_references(references: list) -> list[str]:
        """
        Get only the main author listed inside the reference using regex.

        :parameters
        -----------
        reference - list
            List of the references.

        :returns
        --------
        main_author_references - list[str]
            List of the main authors from the references.

        Format of reference: 
            * authors + year + title + journal edition / book

        Example:
            'Adam, E.J.H. and Adam, S.A. 1994.
            Identification of cytosolic factors required for nuclear location sequence-mediated
            binding to the nuclear envelope. J. Cell Biol. 125:547-555.'

        Try to catch all the sir name + initials:
            * https://regex101.com/r/XhU4Ak/1
            * pattern = r'([ÄÖÜäöüßA-Za-z,\s])+(?:[A-Z]\.){1,10}'
        """
        pattern = r'([ÄÖÜäöüßA-Za-z,\s])+(?:[A-Z]\.){1,10}'
        test_string = 'Adam, E.J.H. and Adam, S.A. 1994. Identification of cytosolic factors required for nuclear location sequence-mediated binding to the nuclear envelope. J. Cell Biol. 125:547-555.'
        # match = re.search(pattern, test_string)
        # print(f"Match: {match}")
        # print(f"Group: {match.group(0)}")
        main_author_references = []
        try:
            for reference in references:
                match = re.search(pattern, reference["Citation"])
                # print(f"Main author: {match.group(0)}")
                main_author_references.append(match.group(0))
        except Exception:
            pass

        return main_author_references


def main():
    """
    Main func.
    """
    data_dir = Path("/data/dataprocessing/NCBI/PubMed/")
    # file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0151.xml"
    # This contains one option for referecing using PMID
    file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0455.xml"

    # This contains the other option for referencing using AUthor ID
    file_2 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # Contains empty reference
    file_3 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # ---- test  1 -----

    # medline_parser = MedlineParser()
    parser = PubmedParser(file_2)

    # records = parser.get_parsed_xml()

    parser.parse_articles()


if __name__ == "__main__":
    main()
