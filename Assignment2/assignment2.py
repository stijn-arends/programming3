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

        parser.add_argument("pubmed_id", action="store", type=str, nargs=1, 
            help="Pubmed ID of the article to harvest for references to download.")

        parser.add_argument("-n", action="store",
                           dest="n", required=False, type=int, default=1,
                           help="Number of peons per client.")

        parser.add_argument("-a", action="store",
                           dest="n", required=False, type=int, default=10,
                           help="Number of references to download concurrently.")

        parser.add_argument('--port_number', "-p", help="The port number that will be used",
                        required=True, type=int)

        parser.add_argument("--host", help="Hosts used, first input is set as the server host", 
            required=True, nargs="+")
        parser.add_argument("-c", help="client mode", action='store_true')
        parser.add_argument("-s", help="server mode", action='store_true')

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

if __name__ == "__main__":
    main()