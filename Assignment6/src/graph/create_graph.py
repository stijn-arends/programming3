"""
Module for creating citation graphs,
either with attributes or without.
"""

# imports
import sys
import os
import time
import pickle
import argparse
from typing import Callable, Any
from functools import wraps
from pathlib import Path
import networkx as nx
import numpy as np

__version__ = "v.01"
__author__ = "Stijn Arends"


def time_func(func) -> Callable:
    """
    A decorator that shows the execution time of another
    function.
    """

    @wraps(func)
    def wrap_func(*args, **kwargs):
        """
        Wrapper function.
        """
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
        print(f'Function {func.__name__!r} executed in {elapsed}\n')
        return result
    return wrap_func


class ArgumentParser:
    """
    Class to parse the input arguments.
    """

    def __init__(self) -> None:
        """Initializer"""
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
            description="Python script creates graphs.",
            epilog="Contact: stijnarend@live.nl")

        parser.version = __version__

        parser.add_argument('--adj_list', action="store",
                    dest="adj_list", required=True,
                    help="File containing an adjacency list.")

        parser.add_argument("--attributes", action="store",
                           dest="attributes", required=True,
                           help="Pickle file containing a dictionary with node attributes.")

        parser.add_argument('--nodes_data', action="store",
                    dest="nodes_data", required=True,
                    help="File containing a list of source nodes to use.")

        parser.add_argument('-o', '--output', action="store",
                    dest="o", required=True,
                    help="Location of output directory.")

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


class Graph:
    """
    Class to create a networkx grapg using an
    adjaceny list and a dictionary with node attribute info
    inisde of a pickle file. 
    """

    def __init__(self, adj_list: Path, attr_dict: Path, node_file: Path) -> None:
        """
        Initializer
        """
        print('Reading adjlist data:')
        self.edge_list: nx.DiGraph = self.read_adjlist_data(adj_list)
        print("Reading attribute info:")
        self.attr_dict: dict = read_pickl(attr_dict)
        print("Reading Node list:")
        self.node_file: list = self.read_nodes_data(node_file)

    def create_graph(self) -> nx.DiGraph:
        """
        Create a citation graph in the form of a networkx
        object (DiGraph) without any attribute information.

        :returns
        --------
        graph - nx.DiGraph
            A DAG citation graph
        """
        graph = self.edge_list
        nodes_list = self.node_file
        graph.add_nodes_from(nodes_list)
        return graph
        
    def create_full_graph(self) -> nx.DiGraph:
        """
        Create a citation graph in the form of a networkx
        object (DiGraph) with node attribute information.

        :returns
        --------
        graph - nx.DiGraph
            A DAG citation graph
        """
        graph = self.edge_list
        nodes_list = self.node_file
        graph.add_nodes_from(nodes_list)
        node_attributes = self.attr_dict
        nx.set_node_attributes(graph, node_attributes)
        return graph

    @staticmethod
    @time_func
    def read_adjlist_data(file) -> nx.DiGraph:
        """
        Read in a adjacency list(i.e., csv file) into a networkx object and
        create a DiGraph.

        :returns
        --------
        graph - nx.DiGraph
            A networkx DiGraph from an adjaceny list.
        """
        graph = nx.read_adjlist(file, create_using=nx.DiGraph())
        return graph

    @staticmethod
    @time_func
    def read_nodes_data(file) -> list[str]:
        """
        Read in a file containing source nodes without
        any information about edges from for example a txt file.

        :returns
        --------
        noddes_list - list[str]
            List of source nodes
        """
        with open(file, 'r') as f:
            nodes_data = f.readlines()

        nodes_list = nodes_data[0].split(' ')
        return nodes_list


@time_func
def to_pickle(data, out_file) -> None:
    """
    Write object to pickle file
    """
    with open(out_file, 'wb') as file_handler:
        pickle.dump(data, file_handler)


def read_pickl(file) -> Any:
    """
    Read in a pickle file.
    """
    with open(file, 'rb') as fh:
        data = pickle.load(fh)

    return data


def main() -> None:
    """main"""
    cla_parser = ArgumentParser()
    adj_list = Path(cla_parser.get_argument('adj_list'))
    attr_dict = Path(cla_parser.get_argument('attributes'))
    nodes_no_refs = Path(cla_parser.get_argument('nodes_data'))
    out_dir = Path(cla_parser.get_argument('o'))
    out_file = out_dir / "citation_graph.pkl"
    out_file_full = out_dir / "citation_graph_full.pkl"

    graph = Graph(adj_list, attr_dict, nodes_no_refs)

    simple_graph = graph.create_graph()
    to_pickle(simple_graph, out_file)

    full_graph = graph.create_full_graph()
    to_pickle(full_graph, out_file_full)
    
    print("Calculate most cited paper:")
    in_degrees = np.array(list(dict(simple_graph.in_degree).values()))
    print(f"Max: {np.max(in_degrees)}")
    index_max = np.argmax(in_degrees)
    print(f"Index of max value: {index_max}")

    print(f"Paper with most citations: {list(simple_graph.in_degree)[index_max]}")


if __name__ == "__main__":
    main()
