import sys
from Bio import Entrez
import multiprocessing as mp
from multiprocessing.managers import BaseManager, SyncManager
import time, queue, os
from pathlib import Path
from typing import Any
from bs4 import BeautifulSoup
import pickle

from arg_parser import ArgumentParser


__author__ = "Stijn Arends"
__version__ = "v0.1"
__data__ = "14-5-2022"

POISONPILL = "MEMENTOMORI"
ERROR = "DOH"
# IP = ''
# PORTNUM = 5381
AUTHKEY = b'whathasitgotinitspocketsesss?'

Entrez.email = "stijnarends@live.nl"
Entrez.api_key = '9f94f8d674e1918a47cfa8afc303838b0408'


def make_server_manager(port, authkey):
    """ Create a manager for the server, listening on the given port.
        Return a manager object with get_job_q and get_result_q methods.
    """
    job_q = queue.Queue()
    result_q = queue.Queue()

    # This is based on the examples in the official docs of multiprocessing.
    # get_{job|result}_q return synchronized proxies for the actual Queue
    # objects.
    class QueueManager(BaseManager):
        pass

    QueueManager.register('get_job_q', callable=lambda: job_q)
    QueueManager.register('get_result_q', callable=lambda: result_q)

    manager = QueueManager(address=('', port), authkey=authkey)
    manager.start()
    print('Server started at port %s' % port)
    return manager


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
        self.get_author_names(paper, pmid)


    def get_author_names(self, paper:bytes, pmid:int) -> None:
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

        # print(f"Authors: {tuple(authors)}")

        output_file = self.out_path / f"{pmid}.authors.pickle"
        self.write_pickle(tuple(authors), output_file)

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
        with open(file, 'wb') as fh:
            pickle.dump(data, fh)


def main():
    # Get passed arguments
    cla_parser = ArgumentParser()
    pmid = cla_parser.get_argument('pubmed_id')
    n_articles= cla_parser.get_argument('a')

    n_peons = cla_parser.get_argument('n')
    port = cla_parser.get_argument('p')
    host = cla_parser.get_argument('host')

    run_server = cla_parser.get_argument('s')
    run_client = cla_parser.get_argument('c')

    print(f"Pubmed ID: {pmid}, articles to download: {n_articles}")
    print(f"Number of peons: {n_peons}")
    print(f"Port number: {port}")
    print(f"host: {host}")
    out_path = Path(os.path.abspath(os.path.dirname(__file__))) / "output"
    download_auths = DownloadAuthorList(n_articles=n_articles, out_path=out_path)

    ref_ids = download_auths.get_id_references(pmid)

    print(f"Reference IDs: {ref_ids}")

    # download_auths.download_paper(ref_ids[0])

    print(f"Server: {run_server}")
    print(f"Client: {run_client}")

    if run_server:
        server = mp.Process(target=runserver, args=(download_auths.download_paper, ref_ids[:n_articles], port))
        server.start()
        time.sleep(1)
        server.join()


    if run_client:
        print('Selected client mode.')
        client = mp.Process(target=runclient, args=(n_peons, host, port))
        client.start()
        client.join()


if __name__ == "__main__":
    main() # 30049270