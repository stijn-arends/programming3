# Assignment 1 of the porgramming 3 course
---

* date: 14-5-2022
* author: Stijn Arends
* assignment [link](https://bioinf.nl/~martijn/master/programming3/assignment1.html)


## Assignment
--- 

A script called "assignment1.py" that given 1 starting article's Pubmed ID, downloads 10 articles referenced by that first article. 
It should do this concurrently from PubMed using the Biopython API.

## Running the script

```bash
python3 assignment1.py 33669327
```

Optionally, the number of articles to download can be specified using the `-n` argument.
```bash
python3 assignment1.py -pmid 33669327 -n 10
```

## Packages
--- 
