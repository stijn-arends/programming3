import sys
from weakref import ref
from Bio import Entrez
from bs4 import BeautifulSoup
import multiprocessing as mp
import argparse
import os

# Assignment: https://bioinf.nl/~martijn/master/programming3/assignment1.html
# Entrenz documentation: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec143

__author__ = "Stijn Arends"
__version__ = "v0.1"
__data__ = "14-5-2022"


class ArticleNotFound(Exception):
    """Exception raised for 404 errors.
    Attributes:
        url -- url of a website
        message -- explanation of the error
    """

    def __init__(self, pmid):
        self.message = f"pmid: {pmid}, is not a valid ID or it does not contain any references."
        super().__init__(self.message)


class DownloadPubmedPapers:
    """
    Download a number of papers that are referenced in a pubmed article. 
    """

    def __init__(self, pmid:str, n_articles:int) -> None:
        self.pmid = pmid
        self.n_articles = n_articles
        Entrez.email = "stijnarends@live.nl"
        Entrez.api_key = '9f94f8d674e1918a47cfa8afc303838b0408'

    def get_id_references(self) -> list:
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
                                id=self.pmid))

        if not results[0]["LinkSetDb"]:
            raise ArticleNotFound(self.pmid)
                                
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
        handle = Entrez.efetch(db="pmc", id=pmid, rettype="XML", retmode="text")
        handle.close()


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
        -v, --version - displays the version of the script
        -h, --help - display the help text

        :returns
        --------
        parser - ArgumentParser
        """
        parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
            description="Python script that downloads 10 papers from the reference section from an article in pubmed.",
            epilog="Contact: stijnarend@live.nl")

        # Set version
        parser.version = __version__

        parser.add_argument("-n", action="store",
                           dest="n", required=False, type=int, default=10,
                           help="Number of references to download concurrently.")

        parser.add_argument('-pmid',
            required=True,
            type=int,
            help='Pubmed ID')

        parser.add_argument('-v',
            '--version',
            help='Displays the version number of the script and exitst',
            action='version')

        return parser

    def get_argument(self, argument_key):
        """
        Method to get an input argument.
        :parameters
        -----------
        argument_key - str 
            Full command line argument (so --config for the configuration file argument).

    
        :returns
        --------
        value - List or boolean
        """
        if self.arguments is not None and argument_key in self.arguments:
            value = getattr(self.arguments, argument_key)
        else:
            value = None
        return value


def main():
    # Get passed arguments
    cla_parser = ArgumentParser()
    pmid = cla_parser.get_argument('pmid')
    n_articles= cla_parser.get_argument('n')

    download_pm = DownloadPubmedPapers(pmid=pmid, n_articles=n_articles)
    refs = download_pm.get_id_references()
    print(refs)

if __name__ == "__main__":
    main()