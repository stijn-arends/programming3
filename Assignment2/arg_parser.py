"""
Module for parsing arguments.
"""

import sys
import argparse
import os
from typing import Any

__author__ = "Stijn Arends"
__version__ = "v0.1"
__data__ = "17-5-2022"


class ArgumentParser:
    """
    Class to parse the input arguments.
    """

    def __init__(self):
        self.parser = self._create_argument_parser()
        # Print help if no arguments are supplied and stop the program
        if len(sys.argv) == 1:
            self.parser.print_help(sys.stderr)
            sys.exit(1)
        self.arguments = self.parser.parse_args()

    @staticmethod
    def _create_argument_parser():
        """
        Create an argument parser.

        :returns
        --------
        parser - ArgumentParser
        """
        parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
            description="Python script that downloads papers from the reference" /
            "section from an article in pubmed.",
            epilog="Contact: stijnarend@live.nl")

        parser.version = __version__

        parser.add_argument("pubmed_id", action="store", type=str, nargs=1,
            help="Pubmed ID of the article to harvest for references to download.")

        parser.add_argument("-n", action="store",
                           dest="n", required=False, type=int, default=1,
                           help="Number of peons per client.")

        parser.add_argument("-a", action="store",
                           dest="a", required=False, type=int, default=10,
                           help="Number of references to download concurrently.")

        parser.add_argument("-p", '--port_number',dest='p',
                        help="The port number that will be used",
                        required=True, type=int)

        parser.add_argument("--host", dest="host",
            help="Hosts used, first input is set as the server host",
            required=True, nargs="?")

        parser.add_argument('-v',
            '--version',
            help='Displays the version number of the script and exitst',
            action='version')

        command_group = parser.add_mutually_exclusive_group()
        command_group.add_argument('-s', dest='s',
            help='Server mode, can\'t be used together with -c',
            action='store_true')

        command_group.add_argument('-c', dest='c',
            help='Client Mode, can\'t be used together with -s',
            action='store_true')

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
