#!/usr/bin/env python3

"""
Module to parse XML files containing the PubMed literature.
"""

# IMPORTS
import datetime
import re
from pathlib import Path
from typing import Any, Tuple

import pandas as pd
from Bio import Entrez

try:
    from pubmed_parser.medline_parser import MedlineParser
except ModuleNotFoundError:
    from medline_parser import MedlineParser

Entrez.email = "stijnarends@live.nl"
Entrez.api_key = "9f94f8d674e1918a47cfa8afc303838b0408"


class PubmedParser:
    """
    A class for parsing pubmed articles in XML format.

    :get_publish_date
    :get_references
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

        return records["PubmedArticle"]

    def get_parsed_xml(self) -> list:
        """
        Get the parsed XML file.

        :returns
        --------
        self.data - list
            Parsed XML data containing info about the articles
        """
        return self.data

    def parse_articles(self) -> pd.DataFrame:
        """
        Create a dictionary of lists and save the results in here.

        :returns
        --------
        df_parsed_article - pd.DataFrame
            The parsed article inside a pandas data frame.
        """
        pmid_list, title_list, author_list, co_author_list = [], [], [], []
        journal_list, key_words_list, language_list, doi_list = [], [], [], []
        pmc_list, publish_date_list, ref_ids_list, ref_title_list = [], [], [], []
        ref_auth_names_list, ref_types = [], []
        for article in self.data:
            medline_parser = MedlineParser(article)
            pmid = medline_parser.get_pmid()
            pmid_list.append(pmid)

            title = medline_parser.get_title()
            title_list.append(title)

            author = medline_parser.get_author()
            author_list.append(author)

            co_authors = medline_parser.get_co_authors()
            co_author_list.append(co_authors)

            journal = medline_parser.get_journal()
            journal_list.append(journal)

            key_words = medline_parser.get_keywords()
            key_words_list.append(key_words)

            language = medline_parser.get_language()
            language_list.append(language)

            doi = medline_parser.get_doi()
            doi_list.append(doi)

            pmc = medline_parser.get_pmc()
            pmc_list.append(pmc)

            publish_date = self.get_publish_date(article)
            publish_date_list.append(publish_date)

            ref_ids, ref_auth_names, ref_titles, ref_type = self.get_references(article)

            ref_ids_list.append(ref_ids)
            ref_title_list.append(ref_titles)
            ref_auth_names_list.append(ref_auth_names)
            ref_types.append(ref_type)

        parsed_data = {
            "pmid": pmid_list,
            "language": language_list,
            "author": author_list,
            "title": title_list,
            "co_authors": co_author_list,
            "journal": journal_list,
            "key_words": key_words_list,
            "doi": doi_list,
            "pmc": pmc_list,
            "publish_date": publish_date_list,
            "ref_ids": ref_ids_list,
            "ref_authors": ref_auth_names_list,
            "ref_titles": ref_title_list,
            "ref_type": ref_types,
        }

        df_parsed_article = pd.DataFrame(parsed_data)
        df_parsed_article["publish_date"] = pd.to_datetime(
            df_parsed_article["publish_date"], errors="coerce"
        )
        df_parsed_article["publish_date"] = df_parsed_article[
            "publish_date"
        ].dt.strftime("%Y-%m-%d")

        return df_parsed_article

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
        publish_date - datetime.date
            The publish_date of publishing, returns empty str when date can't be found.
        """
        try:
            history = article["PubmedData"]["History"]
            pubmed = False
            medline = False
            for i, history_type in enumerate(history):
                if history_type.attributes["PubStatus"] == "pubmed":
                    pubmed = i
                elif history_type.attributes["PubStatus"] == "medline":
                    medline = i

            try:
                if pubmed:
                    year = int(history[pubmed]["Year"])
                    month = int(history[pubmed]["Month"])
                    day = int(history[pubmed]["Day"])
                    publish_date = datetime.date(year, month, day)
                elif medline and not pubmed:
                    year = int(history[medline]["Year"])
                    month = int(history[medline]["Month"])
                    day = int(history[medline]["Day"])
                    publish_date = datetime.date(year, month, day)
                else:
                    publish_date = ""
            except ValueError:  # for incorrect publish_date (i.e., 2000-02-31)
                publish_date = ""
        except (IndexError, KeyError):
            publish_date = ""
        return publish_date

    def get_references(
        self, article: dict
    ) -> Tuple[list[Any], list[list], list[str], str]:
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
        ref_id_list = []
        ref_title_list = []
        ref_author_list = [[""]]
        try:
            reference_list = article["PubmedData"]["ReferenceList"][0]
            references = reference_list["Reference"]
            if len(references[0].keys()) > 1:
                ref_id_list = self._get_pmid_reference(references)
                ref_type = "pmid"
            else:
                ref_author_list, ref_title_list = self._get_authors_title_reference(
                    references
                )
                ref_type = "author"
        except (IndexError, KeyError):
            ref_type = "no-ref"

        return ref_id_list, ref_author_list, ref_title_list, ref_type

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
        pattern = re.compile(r"^\d+$")
        pmid_list = []
        try:
            for reference in references:
                article_id_list = reference["ArticleIdList"][0]
                pmid = (
                    article_id_list
                    if article_id_list.attributes["IdType"] == "pubmed"
                    else ""
                )
                match = pattern.match(pmid)  # check if pmid only consits of digits
                if match:
                    pmid_list.append(int(pmid))
        except (IndexError, KeyError):
            pass
        return pmid_list

    def _get_authors_title_reference(
        self, references: list
    ) -> Tuple[list[list[str]], list[str]]:
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
        pattern_authors = r"([ÄÖÜäöüßA-Za-z,\s.])+(?:[A-Z]\.|[A-Z],){1,10}"
        pattern_et_al = r"et al."
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
        names = names.replace("and ", "").replace(".", "").replace(",", "").rstrip()
        split_names = names.split(" ")
        correct_names = []
        for i in range(1, len(split_names), 2):
            corr_initals = ".".join(split_names[i]) + "."
            correct_names.append(split_names[i - 1] + ", " + corr_initals)
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
        pattern_year = (
            r"^[\d\.|\d\.\(\)\s]+"  # capture possible year that is infront title
        )
        pattern_title = r"^.*?(?=\.)"  # Grab the first sentence until a dot

        # Remove year infront of title if it excists
        final_str = re.sub(pattern_year, "", reference)

        # Grab title from sentence
        title_match = re.match(pattern_title, final_str)
        if title_match:
            return title_match.group(0)

        return ""


def main():
    """
    Main func.
    """
    # data_dir = Path("/data/dataprocessing/NCBI/PubMed/")

    # This contains one option for referecing using PMID
    # file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0455.xml"

    # This contains the other option for referencing using AUthor ID
    # file_2 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # # Contains empty reference
    # file_3 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # # File that contains a date error - problem solved
    # file_4 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0379.xml"

    # File that contains an error whilst trying to get the pmid reference - problem solved
    file_6 = Path("/data/dataprocessing/NCBI/PubMed/pubmed21n1049.xml")

    # ---- test  1 -----

    parser = PubmedParser(file_6)

    # records = parser.get_parsed_xml()

    # df = parser.parse_articles()

    # print(df.head(10))

    # print(df.ref_titles)

    # title_refs = df[df["ref_titles"].str.len() != 0]
    # print(title_refs)

    spark_df = parser.parse_articles()
    print(spark_df.printSchema())
    print(spark_df.show())

    author_refs = spark_df.filter(spark_df.ref_type == "author")
    print(author_refs.show())

    # test = spark_df.select('ref_titles').take(5)
    # print(type(test[0][0]))
    # print(test)

    # example_list = title_refs.iloc[0, -2]
    # print(f"Example list: {example_list}")
    # print(type(example_list))

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
