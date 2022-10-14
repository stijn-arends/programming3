#!/usr/bin/env python3

"""
Module to parse medline data.
"""

#IMPORTS
# import glob
import re
from typing import Any

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
        except (IndexError, KeyError):
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
        except (IndexError, KeyError):
            pmid = ""
        return pmid

    @staticmethod
    def _search_pmid_article_list(article_id_list: list) -> str:
        """
        Search for the PMID inside of the ArticleIdList which is a key inside
        of the PubmedData.

        :parameter
        ----------
        article_id_list - list
            List of various IDs (i.e., pmid, PMC)

        :returns
        --------
        pmid - str
            PMID of an article.
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
        except (IndexError, KeyError):
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
                    if "LastName" in author and "Initials" in author:
                        last_name = author["LastName"]
                        initials = author["Initials"]
                        # add . between and after the initals
                        if '.' not in initials:
                            initials = ".".join(initials) + '.'
                        co_author_names.append(str(last_name + ", " + initials))
        except (IndexError, KeyError):
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
        except (IndexError, KeyError):
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
        except (IndexError, KeyError):
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
        except (IndexError, KeyError):
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
        except (KeyError, IndexError):
            doi = ""
        return doi

    @staticmethod
    def _search_doi(doi_list: list) -> str:
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
        except (IndexError, KeyError):
            pmc = ""
        return pmc