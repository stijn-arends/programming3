from Bio import Entrez
from bs4 import BeautifulSoup
import multiprocessing as mp

# Assignment: https://bioinf.nl/~martijn/master/programming2/assignment1.html
# Entrenz documentation: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec143


Entrez.email = "stijnarends@live.nl"  # Always tell NCBI who you are
handle = Entrez.einfo()
print(handle)
#result = handle.read()
# handle.close()

 
# Passing the stored data inside
# the beautifulsoup parser, storing
# the returned object
# Bs_data = BeautifulSoup(result, "xml")
 
# Finding all instances of tag
# `unique`
# b_unique = Bs_data.find_all('DbName')

 
# print(b_unique)


# ----
# Use entrenz.efetch(pmid) to get the paper in xml format.
# ----


def retrieve_ref_ids(pmid):
    ids = []
    record = Entrez.read(Entrez.elink(dbfrom="pubmed", id=pmid))


    # Get the references of the paper
    index = 5
    for i, rec in enumerate(record[0]["LinkSetDb"]):
        # print(rec)
        if  rec['LinkName'] == "pubmed_pubmed_refs":
            index = i
    
    for link in record[0]["LinkSetDb"][index]["Link"]:
        ids.append(link["Id"])

    return ids


def fetch_paper(pmid):
    paper = Entrez.efetch(db="pubmed", id=pmid, retmode='xml').read() # retmode='txt', rettype='xml'
    write_xml(paper, pmid)


def write_xml(data, name):
    bs_data = BeautifulSoup(data, "xml")

    f = open(f'{name}.xml', "w", encoding="utf-8")
    f.write(bs_data.prettify())
    f.close()


def main():
    # pmid = "35223391" # "19304878" # "33669327" # 
    pmid = "19304878"
    # results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", id=pmid))
    # pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
    # print(pmc_ids)
    # print(len(pmc_ids))


    ids = retrieve_ref_ids(pmid)

    print(ids)
    print(len(ids))

    test_id = ids[1]
    print(test_id)

    # paper = Entrez.efetch(db="pubmed", id=test_id, retmode='xml').read()

    # bs_data = BeautifulSoup(paper, "xml")

    # print(str(bs_data))

    # f = open(f'{test_id}.xml', "w", encoding="utf-8")
    # f.write(bs_data.prettify())
    # f.close()

    # cpus = mp.cpu_count()

    # with mp.Pool(4) as pool:
    #     pool.map(fetch_paper, ids[:10])
    for id in ids[:10]:
        fetch_paper(id)

if __name__ == "__main__":
    main()