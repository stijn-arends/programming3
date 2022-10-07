#!/usr/bin/env python3

"""
Module to parse XML files containing the PubMed literature.
"""

#IMPORTS
from bs4 import BeautifulSoup
from Bio import Entrez, Medline


Entrez.email = "stijnarends@live.nl"
Entrez.api_key = '9f94f8d674e1918a47cfa8afc303838b0408'

def read_xml(file):
    with open(file, "br") as handle:
        records = Entrez.read(handle)
    return records['PubmedArticle']

def read_bs4(file):
    with open(file, 'br') as f:
        file_contents = f.read() 

    # 'xml' is the parser used. For html files, which BeautifulSoup is typically used for, it would be 'html.parser'.
    soup = BeautifulSoup(file_contents, 'xml')
    print(soup)


def main():
    data_dir = "/data/dataprocessing/NCBI/PubMed/"
    # file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0151.xml"
    file = "/data/dataprocessing/NCBI/PubMed/pubmed21n0455.xml"

    records = read_xml(file)
    print(len(records))
    print(records[0].keys())
    print(records[0]["MedlineCitation"].keys())
    print(records[0]['PubmedData'].keys())

    print(f"Reference list: {records[0]['PubmedData']['ReferenceList']}")
    print(f"Title: {records[0]['MedlineCitation']['Article']['ArticleTitle']}")
    print(len(records[0]))

    print("\n")

    #read_bs4(file)
    print(records[1]["MedlineCitation"].keys())
    print(records[1]['PubmedData'].keys())

    print(f"Reference list: {records[100]['PubmedData']['ReferenceList']}")
    print(f"Title: {records[100]['MedlineCitation']['Article']['ArticleTitle']}")
    print(len(records[100]))

    count = 0
    for article in records:
        if article['PubmedData']['ReferenceList'] and count < 3:
            print(f"\nNumber: {count + 1}\n")
            # print(article['PubmedData']['ReferenceList'])
            print(article['PubmedData']['ReferenceList'])
            print(type(article['PubmedData']['ReferenceList']))
            print()
            print(article['PubmedData']['ReferenceList'][0].keys())
            # print(len(article['PubmedData']['ReferenceList']))
            count += 1

    count2 = 0
    count3 = 0
    for article in records:
        if article['PubmedData']['ReferenceList']:
            if len(article['PubmedData']['ReferenceList'][0]["Reference"][0]) == 2:
                count2 += 1
            else:
                count3 += 1

    print(f"Total number of articles: {len(records)}")
    print(f"Number of articles with 2 keys: {count2}")
    print(f"Other: {count3}")


if __name__ == "__main__":
    main()