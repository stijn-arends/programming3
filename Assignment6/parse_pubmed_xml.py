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
from typing import Any, Tuple
# from bs4 import BeautifulSoup
from Bio import Entrez
import pandas as pd
import ast


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
        pmid_list, title_list, author_list, co_author_list = [], [], [], []
        journal_list, key_words_list, language_list, doi_list = [], [], [], []
        pmc_list, publish_date_list, ref_ids_list, ref_title_list = [], [], [], []
        for article in self.data: # [12:13]
            medline_parser = MedlineParser(article)
            pmid = medline_parser.get_pmid()
            pmid_list.append(pmid)
            # print(f"\nPubmed ID: {pmid}")
            title = medline_parser.get_title()
            title_list.append(title)
            # print(f"Title: {title}")
            author = medline_parser.get_author()
            author_list.append(author)
            # print(f"Author: {author}")
            co_authors = medline_parser.get_co_authors()
            co_author_list.append(co_authors)
            # print(f"Co-authors: {co_authors}")
            journal = medline_parser.get_journal()
            journal_list.append(journal)
            # print(f"Journal: {journal}")
            key_words = medline_parser.get_keywords()
            key_words_list.append(key_words)
            # print(f"Key words: {key_words}")
            language = medline_parser.get_language()
            language_list.append(language)
            # print(f"Language: {language}")
            doi = medline_parser.get_doi()
            doi_list.append(doi)
            # print(f"DOI: {doi}")
            pmc = medline_parser.get_pmc()
            pmc_list.append(pmc)
            # print(f"PMC: {pmc}")
            publish_date = self.get_publish_date(article)
            publish_date_list.append(publish_date)
            # print(f"Publish date: {publish_date}")
            ref_ids, ref_titles = self.get_references(article)
            # print(f"Reference ID list: {ref_ids}")
            ref_ids_list.append(ref_ids)
            ref_title_list.append(ref_titles)

        parsed_data = {"pmid": pmid_list, "language": language_list,
                    "author": author_list, "title": title_list,
                    "co_authors": co_author_list, "journal": journal_list,
                    "key_words": key_words_list, "doi": doi_list, "pmc": pmc_list,
                    "publish_date": publish_date_list, "ref_ids": ref_ids_list,
                    "ref_titles": ref_title_list}
        # Perhaps write the df out as pickle files
        df = pd.DataFrame(parsed_data)
        # print(df.head(10))
        # print(df.shape)
        # print(sum(df.author == ""))
        df["publish_date"] = df["publish_date"].astype("datetime64")
        return df


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
            history = article['PubmedData']['History']
            year = int(history[0]['Year'])
            month = int(history[0]['Month'])
            day = int(history[0]['Day'])
            # print(f"{article['PubmedData']['History']}")
            # print(f"Year: {year}, month: {month}, day: {day}\n")
            try:
                date = datetime.date(year, month, day)
            except ValueError:
                year = int(history[1]['Year'])
                month = int(history[1]['Month'])
                day = int(history[1]['Day'])
                date = datetime.date(year, month, day)
        except (IndexError, KeyError):
            date = ""
        return date

    def get_references(self, article: dict) -> Tuple[list[Any], list[str]]:
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
        ref_title_list - list[str]
            List of the titles of the references.

        :see also
        ---------
        PubmedParser._get_pmid_reference()
        PubmedParser._get_author_title_reference()
        PubmedParser._get_author_names()
        PubmedParser._get_title_reference()
        """
        try:
            reference_list = article['PubmedData']['ReferenceList'][0]
            references = reference_list['Reference']
            if len(references[0].keys()) > 1:
                ref_id_list = self._get_pmid_reference(references)
                ref_title_list = []
            else:
                ref_id_list, ref_title_list = self._get_authors_title_reference(references)
        except (IndexError, KeyError):
            ref_id_list = []
            ref_title_list = []

        return ref_id_list, ref_title_list

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
        pmid_list - list[int]
            List of the PMID from the references.
        """
        pattern = re.compile(r'^\d+$')
        pmid_list = []
        try:
            for reference in references:
                article_id_list = reference["ArticleIdList"][0]
                pmid = article_id_list if article_id_list.attributes["IdType"] == "pubmed" else ""
                match = pattern.match(pmid) # check if pmid only consits of digits
                if match:
                    pmid_list.append(int(pmid))
        except (IndexError, KeyError):
            pass
        return pmid_list

    def _get_authors_title_reference(self, references: list) -> Tuple[list[str], list[str]]:
        """
        Get the names of the authors and the title for each reference using regex.
        First, remove 'et al' from the reference. Next, extract all authors from the reference.
        If there was a match proceed to process the string of authors to seperate them into a list.
        Remove the authors from the original reference, then, grab the title from the reference.

        regex test; extract author names: https://regex101.com/r/X3ov1K/1
        regex test; match et al: https://regex101.com/r/R7HDoe/1

        Format of reference:
            * authors + year + title + journal edition / book

        Example:
            'Adam, E.J.H. and Adam, S.A. 1994.
            Identification of cytosolic factors required for nuclear location sequence-mediated
            binding to the nuclear envelope. J. Cell Biol. 125:547-555.'

        :parameters
        -----------
        reference - list
            List of the references.

        :returns
        --------
        authors_references - list[str]
            List of the authors from the references.
        titles_references - list[str]
            List of the titles from the references.

        :see also
        ---------
        PubmedParser._get_author_names()
        PubmedParser._get_titile_reference()
        """
        pattern_authors = r"([ÄÖÜäöüßA-Za-z,\s.])+(?:[A-Z]\.|[A-Z],){1,10}" # r"([ÄÖÜäöüßA-Za-z,\s.])+(?:[A-Z]\.){1,10}"
        pattern_et_al = r'et al.'
        authors_references = []
        titles_references = []

        try: 
            for reference in references:
                citation = reference["Citation"]
                citation = re.sub(pattern_et_al, "", citation)
                match = re.match(pattern_authors, citation)
                if match:
                    names = match.group(0)

                    corr_names = self._get_author_names(names)
                    authors_references.append(corr_names)

                    # Remove the author names from the reference
                    new_str = re.sub(pattern_authors, "", citation)
                    
                    # Get the title
                    title = self._get_title_reference(new_str)
                    titles_references.append(title)
                else:
                    title = self._get_title_reference(citation)
                    authors_references.append([])
                    titles_references.append(title)
        except (IndexError, KeyError):
            pass

        return authors_references, titles_references

    @staticmethod
    def _get_author_names(names: str) -> list[str]:
        """
        Parse a string of multiple author to seperate them into
        a list.

        :parameter
        ----------
        names - str
            Names of authors; sirname + initials

        :returns
        --------
        correct_names - list[str]
            List of authors sirname and initals.
        """
        names = names.replace('and ', '').replace('.', '').replace(',', '').rstrip()
        names = names.split(' ')
        correct_names = []
        for i in range(1, len(names), 2):
            corr_initals = ".".join(names[i]) + '.'
            correct_names.append(names[i -1] + ", " + corr_initals)
        return correct_names

    @staticmethod
    def _get_title_reference(reference: str) -> str:
        """
        Extract the title from a reference.


        regex test: extract year from beginning of sentence: https://regex101.com/r/k5VMfy/1
        regex test: match everythin until first dot: https://regex101.com/r/iw8o2y/1

        :parameter
        ----------
        reference - str
            A reference

        :returns
        --------
        title - str
            The extracted title.
        """
        pattern_year = r'^[\d\.|\d\.\(\)\s]+' # capture possible year that is infront title
        pattern_title = r'^.*?(?=\.)' # Grab the first sentence until a dot

        # Remove year infront of title if it excists
        final_str = re.sub(pattern_year, "", reference)

        # Grab title from sentence
        title_match = re.match(pattern_title, final_str)
        if title_match:
            return title_match.group(0)
        else:
            return ""


def main():
    """
    Main func.

    Something goes wrong with the date with file pubmed21n0379.xml:

    When writing out the CSV file, everything becomes a str.
    We can't have this because we would like to keep the lists inside the
    dataframe. This stackoverflow post discusses this:
        * https://stackoverflow.com/questions/48250995/write-lists-to-pandas-dataframe-to-csv-read-dataframe-from-csv-and-convert-to-l
    """
    data_dir = Path("/data/dataprocessing/NCBI/PubMed/")

    # This contains one option for referecing using PMID
    # file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0455.xml"

    # This contains the other option for referencing using AUthor ID
    file_2 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # Contains empty reference
    file_3 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # File that contains a date error - problem solved
    file_4 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0379.xml"

    # File that contains an error that occurs in _get_main_author_refences - no longer using this function
    file_5 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0632.xml"

    # File that contains an error whilst trying to get the pmid reference - problem solved
    file_6 = "/data/dataprocessing/NCBI/PubMed/pubmed21n1049.xml"

    # ---- test  1 -----

    parser = PubmedParser(file_6)

    # records = parser.get_parsed_xml()

    df = parser.parse_articles()

    print(df.head(10))

    print(df.ref_titles)

    title_refs = df[df["ref_titles"].str.len() != 0]
    print(title_refs)

    example_list = title_refs.iloc[0, -2]
    print(f"Example list: {example_list}")
    print(type(example_list))

    # example_list_title = title_refs.iloc[0, -1]
    # print(f"Example list: {example_list_title}")
    # print(type(example_list_title))

    # ----- test all files ----

    # out_dir = Path("/commons/dsls/dsph/2022/parsed_pubmed_articles")
    # out_dir.mkdir(exist_ok=True)

    # for i, file in enumerate(data_dir.glob("*.xml")):
    #     # if i > 1047:
    #     print(f"Parsing file: n: {i+1}, name: {file.stem}")
    #     parser = PubmedParser(file)
    #     df = parser.parse_articles()
    #     out_file = out_dir / (file.stem + ".csv")
    #     df.to_csv(out_file, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()
