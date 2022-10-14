"""
Module for parsing arguments.
"""

import sys
import argparse
import os
from pathlib import Path
from typing import Any
from typing_extensions import Required

__author__ = "Stijn Arends"
__version__ = "v0.1"
__data__ = "14-10-2022"


class ArgumentParser:
    """
    Class to parse the input arguments.
    """

    def __init__(self) -> None:
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
            description="Python script that parses pubmed articles",
            epilog="Contact: stijnarend@live.nl")

        parser.version = __version__

        parser.add_argument("-d", '--data_dir', action="store",
                    dest="d", required=True, type=str,
                    help="Location of the pubmed data.")

        parser.add_argument("-n", action="store",
                           dest="n", required=False, type=int, default=1,
                           help="Number of peons per client.")

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

        command_group = parser.add_mutually_exclusive_group(required=True)
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

class CLIArgValidator:
    """
    Class to check if arguments are valid.
    """

    def validate_input_file(self, input_path: str) -> None:
        """
        Validate the input files by checking if they actually exists.
        :parameters
        -----------
        input_path - str
            Path to a file
        """
        input_path = Path(input_path)
        self._validate_input_exists(input_path)
        self._check_xml_files(input_path)

    @staticmethod
    def _validate_input_exists(input_path: Path) -> None:
        """
        Check if a file exists.

        :parameters
        -----------
        input_path - Path
            Path to a file
        """
        if not input_path.exists():
            raise FileNotFoundError(f'Input file does not exist!: {input_path}'\
                '\nPlease check if directory was specified correctly.')

    @staticmethod
    def _check_xml_files(input_path: Path) -> None:
        """
        Check if the specified data dir contains any XML files.

        :parameters
        -----------
        input_path - Path
            Path to a file
        """
        xml_files_count = 0
        for file in input_path.iterdir():
            if file.stem.endswith("xml"):
                xml_files_count += 1

        if xml_files_count == 0:
            raise ValueError(f"No XML files were found inside: {input_path}")
