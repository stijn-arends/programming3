#!/usr/bin/env python3
"""
From a given pubmed IDs retrieve the pubmed IDs of the references
from that article.

A server puts those IDs in a job queue, a client will connect to
the server and process an ID by downloading the authors from the article
and saving it in a pickle file.
"""

import multiprocessing as mp
import time
from pathlib import Path
from typing import Any
import pickle
from Bio import Entrez

from arg_parser import ArgumentParser
from server_side import ServerSide
from client_side import ClientSide


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


class DownloadPubmedInfo:
    """
    Download a number of papers that are referenced in a pubmed article.
    """

    def __init__(self, n_articles:int, out_path: Path) -> None:
        self.n_articles = n_articles
        self.out_path = out_path
        self.make_data_dir(out_path)

    def make_data_dir(self, path: Path) -> None:
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
            print(f"[{self.make_data_dir.__name__}] {path} already exists.")

    @staticmethod
    def get_id_references(pmid) -> list:
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

    def get_info(self, pmid:int) -> None:
        """
        Get pubmed info.

        :parameters
        ----------
        pmid - int
            Pubmed ID
        """
        self.get_authors(pmid)
        self.download_paper(pmid)

    def get_authors(self, pmid:int) -> None:
        """
        Get the names of the authors of a paper.

        :parameter
        ----------
        pmid - int
            Pubmed ID
        """
        try:
            with Entrez.esummary(db="pubmed", id=pmid) as handle:
                info = Entrez.read(handle)
                authors = tuple(info[0]['AuthorList'])
        except RuntimeError:
            print(f"Nothing found for pumed ID: {pmid}")
            authors = (None)

        output_file = self.out_path / f"{pmid}.authors.pickle"
        self.write_pickle(authors, output_file)

    @staticmethod
    def write_pickle(data: Any, file: Path) -> None:
        """
        Write data to a pickle file

        :parameters
        -----------
        data - Any
            Any data that can be pickled
        file - Path
            Name and location of the file to write to
        """
        with open(file, 'wb') as file_handle:
            pickle.dump(data, file_handle)

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
        self.write_out_paper(paper, pmid)

    def write_out_paper(self, data:str, pmid:str) -> None:
        """
        Write out a paper in XML format.

        :parameters
        -----------
        data - str
            Content of a paper in a binary string.
        pmid - int
            Pubmed ID
        """
        with open(self.out_path / f'{pmid}.xml', 'wb') as file:
            file.write(data)


def main():
    """
    Get the supplied arguments and run the server and client code.
    """
    start_time = time.time()
    # Get passed arguments
    cla_parser = ArgumentParser()
    pmid = cla_parser.get_argument('pubmed_id')
    n_articles= cla_parser.get_argument('a')

    n_peons = cla_parser.get_argument('n')
    port = cla_parser.get_argument('p')
    host = cla_parser.get_argument('host')

    server_mode = cla_parser.get_argument('s')
    client_mode = cla_parser.get_argument('c')

    download_auths = DownloadPubmedInfo(n_articles=n_articles,
        out_path=Path(__file__).parent.absolute() / 'output')

    ref_ids = download_auths.get_id_references(pmid)

    if n_articles > len(ref_ids):
        print(f"There are only {len(ref_ids)} articles - specified: {n_articles}")
        n_articles = len(ref_ids)

    if server_mode:
        server_side = ServerSide(ip_adress=host, port=port,
            auth_key=b'whathasitgotinitspocketsesss?',
            poison_pill="MEMENTOMORI")
        server = mp.Process(target=server_side.run_server, args=(download_auths.get_info,
            ref_ids[:n_articles]))
        server.start()
        time.sleep(1)
        server.join()

    if client_mode:
        print('Selected client mode.')
        client_side = ClientSide(ip_adress=host, port=port,
            auth_key=b'whathasitgotinitspocketsesss?',
            poison_pill="MEMENTOMORI")
        client = mp.Process(target=client_side.run_client, args=(n_peons,))
        client.start()
        client.join()

    print(f"\n--- {time.time() - start_time} seconds ---")


if __name__ == "__main__":
    main() # 30049270
