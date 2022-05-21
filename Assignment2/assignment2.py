from Bio import Entrez
import multiprocessing as mp
import time, os
from pathlib import Path
from typing import Any
import pickle

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

    server_mode = cla_parser.get_argument('s')
    client_mode = cla_parser.get_argument('c')

    print(f"Pubmed ID: {pmid}, articles to download: {n_articles}")
    print(f"Number of peons: {n_peons}")
    print(f"Port number: {port}")
    print(f"host: {host}")
    
    out_path = Path(os.path.abspath(os.path.dirname(__file__))) / "output"
    download_auths = DownloadAuthorList(n_articles=n_articles, out_path=out_path)

    ref_ids = download_auths.get_id_references(pmid)

    print(f"Reference IDs: {ref_ids}")


    auth_key= b'whathasitgotinitspocketsesss?'
    poison_pill = "MEMENTOMORI"

    if server_mode:
        server_side = ServerSide(ip=host, port=port, auth_key=auth_key, poison_pill=poison_pill)
        server = mp.Process(target=server_side.run_server, args=(download_auths.get_authors, ref_ids[:n_articles]))
        server.start()
        time.sleep(1)
        server.join()


    if client_mode:
        print('Selected client mode.')
        client_side = ClientSide(ip=host, port=port, auth_key=auth_key, poison_pill=poison_pill)
        client = mp.Process(target=client_side.run_client, args=(n_peons,))
        client.start()
        client.join()


if __name__ == "__main__":
    main() # 30049270