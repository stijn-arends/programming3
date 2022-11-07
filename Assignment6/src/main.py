"""
The main module.
"""

__author__ = "Stijn Arends"
__version__ = "v.01"

import json
import multiprocessing as mp
import time
# IMPORTS
from pathlib import Path

from arg_parser import ArgumentParser, CLIArgValidator
from parallel_computing.client_side import ClientSide
from parallel_computing.server_side import ServerSide
from pubmed_parser.parse_pubmed_xml import PubmedParser


def parse_file(file: Path, out_dir: Path) -> None:
    """
    Parse a pubmed article and write the contents to a CSV file.

    :parameters
    -----------
    file - Path
        pubmed file in xml format
    out_dir - Path
        output directory to store the results
    """
    parser = PubmedParser(file)
    data_frame = parser.parse_articles()
    out_file = out_dir / (file.stem + ".json")
    json_data = json.loads(data_frame.to_json(orient="records", date_format="iso"))

    with open(out_file, "w", encoding="utf-8") as file_handler:
        file_handler.write(json.dumps(json_data))


def make_output_dir(path: Path) -> None:
    """
    Create a directory (if it does not exsit yet) to store the
    data.

    :parameter
    ----------
    path - Path
        Path to a directory

    :Excepts
    --------
    FileExistsError
        The directory already exists
    """
    try:
        path.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(f"[{make_output_dir.__name__}] {path} already exists.")


def main() -> None:
    """
    Main script.
    """
    start_time = time.time()
    cla_parser = ArgumentParser()
    cla_validator = CLIArgValidator()
    data_dir = cla_parser.get_argument("d")

    if data_dir:
        data_dir = Path(data_dir)
        data_files = list(data_dir.glob("*.xml"))
        cla_validator.validate_input_file(data_dir)

    output_dir = cla_parser.get_argument("o")
    if output_dir:
        output_dir = Path(output_dir)
        make_output_dir(output_dir)

    n_peons = cla_parser.get_argument("n")
    port = cla_parser.get_argument("p")
    host = cla_parser.get_argument("host")

    server_mode = cla_parser.get_argument("s")
    client_mode = cla_parser.get_argument("c")

    # ---- Test server client mode ------

    if server_mode:
        server_side = ServerSide(
            ip_adress=host,
            port=port,
            auth_key=b"whathasitgotinitspocketsesss?",
            poison_pill="MEMENTOMORI",
        )
        server = mp.Process(
            target=server_side.run_server, args=(parse_file, data_files, output_dir)
        )
        server.start()
        time.sleep(1)
        server.join()

    if client_mode:
        print("Selected client mode.")
        client_side = ClientSide(
            ip_adress=host,
            port=port,
            auth_key=b"whathasitgotinitspocketsesss?",
            poison_pill="MEMENTOMORI",
        )
        client = mp.Process(target=client_side.run_client, args=(n_peons,))
        client.start()
        client.join()

    elapsed_time = time.time() - start_time
    days = 0
    if elapsed_time >= 86400:
        days = int(elapsed_time / 86400)
    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print(f"\n---- Time: {days}:{elapsed} ----")
    # ----- End server client mode ------


if __name__ == "__main__":
    main()
