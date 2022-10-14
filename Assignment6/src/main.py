"""
The main module.
"""

__author__ = "Stijn Arends"
__version__ = "v.01"

# IMPORTS
from pathlib import Path

from pubmed_parser.parse_pubmed_xml import PubmedParser

def main():
    data_dir = Path("/data/dataprocessing/NCBI/PubMed/")

    # This contains one option for referecing using PMID
    # file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0455.xml"

    # This contains the other option for referencing using AUthor ID
    file_2 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # Contains empty reference
    file_3 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0591.xml"

    # File that contains a date error - problem solved
    file_4 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0379.xml"

    # File that contains an error that occurs in _get_main_author_refences - no longer using this function
    file_5 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0632.xml"

    # File that contains an error whilst trying to get the pmid reference - problem solved
    file_6 = "/data/dataprocessing/NCBI/PubMed/pubmed21n1049.xml"

    # Another datetime error in both files - year is equal to 1 which causes an outofdatebound error.
    file_7 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0542.xml"
    file_8 = "/data/dataprocessing/NCBI/PubMed/pubmed21n0943.xml"
    # pd.Timestamp('15000101')

    # ---- test  1 -----

    # parser = PubmedParser(file_8)

    # # records = parser.get_parsed_xml()

    # df = parser.parse_articles()

    # print(df.head(10))

    # # print(df.ref_titles)

    # title_refs = df[df["ref_titles"].str.len() != 0]
    # print(title_refs)

    # example_list = title_refs.iloc[0, -2]
    # print(f"Example list: {example_list}")
    # print(type(example_list))

    # ----- test all files ----

    out_dir = Path("/commons/dsls/dsph/2022/parsed_pubmed_articles")
    out_dir.mkdir(exist_ok=True)

    for i, file in enumerate(data_dir.glob("*.xml")):
        if i > 941:
            print(f"Parsing file: n: {i+1}, name: {file.stem}")
            parser = PubmedParser(file)
            df = parser.parse_articles()
            out_file = out_dir / (file.stem + ".csv")
            df.to_csv(out_file, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()