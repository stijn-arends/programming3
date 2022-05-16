import sys
from Bio import Entrez
import multiprocessing as mp
import argparse
import os
from pathlib import Path
from typing import Any

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


class DownloadPubmedPapers:
    """
    Download a number of papers that are referenced in a pubmed article. 
    """

    def __init__(self, n_articles:int, out_path) -> None:
        self.n_articles = n_articles
        self.out_path = out_path
        self.make_data_dir(self.out_path)

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



class ArgumentParser:
    """
    Class to parse the input arguments.
    """

    def __init__(self):
        parser = self._create_argument_parser()
        # Print help if no arguments are supplied and stop the program
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)
        self.arguments = parser.parse_args()

    @staticmethod
    def _create_argument_parser():
        """
        Create an argument parser.
        
        :arguments
        ----------
        -pmid - pubmed id
        -n - number of articles to download. 
        -v, --version - displays the version of the script
        -h, --help - display the help text

        :returns
        --------
        parser - ArgumentParser
        """
        parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
            description="Python script that downloads papers from the reference section from an article in pubmed.",
            epilog="Contact: stijnarend@live.nl")

        parser.version = __version__

        parser.add_argument("-n", action="store",
                           dest="n", required=False, type=int, default=10,
                           help="Number of references to download concurrently.")

        parser.add_argument("pubmed_id", action="store", type=str, nargs=1, help="Pubmed ID of the article to harvest for references to download.")

        parser.add_argument('-v',
            '--version',
            help='Displays the version number of the script and exitst',
            action='version')

        return parser

    def get_argument(self, argument_key: str) -> Any:
        """
        Method to get an input argument.

        :parameters
        -----------
        argument_key - str
            Name of command line argument.

        :returns
        --------
        value - Any
            Value of a command line argument
        """
        if self.arguments is not None and argument_key in self.arguments:
            value = getattr(self.arguments, argument_key)
        else:
            value = None
        return value


def main():
    # Get passed arguments
    cla_parser = ArgumentParser()
    pmid = cla_parser.get_argument('pubmed_id')
    n_articles= cla_parser.get_argument('n')

    out_path = Path(os.path.abspath(os.path.dirname(__file__))) / "output"

    download_pm = DownloadPubmedPapers(n_articles=n_articles, out_path=out_path)
    ref_ids = download_pm.get_id_references(pmid=pmid)

    if n_articles > len(ref_ids):
        print(f"There are only {len(ref_ids)} articles - specified: {n_articles}")
        n_articles = len(ref_ids)

    with mp.Pool(mp.cpu_count()) as p:
        p.map(download_pm.download_paper, ref_ids[0:n_articles])

    print(f"Success! Downloaded {n_articles} articels.")

if __name__ == "__main__":
    main()