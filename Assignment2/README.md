# Assignment 2 of the porgramming 3 course
* * *

* Date: 21-5-2022
* Author: Stijn Arends
* Assignment: [link](https://bioinf.nl/~martijn/master/programming3/assignment2.html)

* * *
## Assignment

Your script needs to download "numberofarticles" which are referenced by "STARTINGPUBMEDID" using "numberofchildren" processes on the "hosts" computers. 
(So, like assignment1, but now from multiple computers). Each child can write to the same folder on the NFS filesystem (your home, or somewhere on the /commons).

Your script needs to analyze the XML of each of the references further to extract all the authors of the article. It should save the authors in a 
Python tuple and use the Pickle module to save it to the disk as `output/PUBMED_ID.authors.pickle` where PUBMEDID is of course the pubmed ID of the article in question.

## Running the script

```bash
python3 assignment2.py -n <number_of_peons_per_client> [-c | -s] --port <portnumber> --host <serverhost> -a <number_of_articles_to_download> STARTING_PUBMED_ID
```

* * *
## Packages

| Package                                                           | Version        |
| ----------------------------------------------------------------- | :------------: |
| [biopython](https://biopython.org/)                               | `1.79`         |
| [pathlib](https://pathlib.readthedocs.io/en/0.5/l)                | `1.0.1`        |
