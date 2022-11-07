<div id="top"></div>

<!-- [![Contributors][contributors-shield]][contributors-url] -->
<!-- [![MIT License][license-shield]][license-url] -->

<div align="center">
<h1 align="center">Investigation of the scientific literature</h1>
<h2 align="center">Assignment 6 - programming 3</h2>
  <a href="https://pubmed.ncbi.nlm.nih.gov/">
    <img src="figures/PubMed-Logo.png" alt="Logo" width="30%" height="90">
  </a>
  <a href="https://spark.apache.org/docs/latest/api/python/index.html">
    <img src="figures/Apache_Spark_logo.png" alt="spark" width="30%" height="90" align='right'>
  </a>
   <a href="https://networkx.org/">
    <img src="figures/networkx_logo.svg" alt="spark" width="30%" height="90" align='left'>
  </a>
</div>
<br>

* * *
# About this project
This is the [final assignment](https://bioinf.nl/~martijn/master/programming3/assignment6.html) of a series of in total 6 assignment that were assigned for the course programming 3 which is all about working with big data.

The goal of this project is to investigate the scientific literature. This will be done by parsing and analyzing the entire PubMed literature database in XML format. Graph theory and parallel processing need to be applied in order to process and analyze this data. 

Furthermore, we wish to answer these specific questions:
1. How large a group of co-authors does the average publication have?
2. Do authors mostly publish using always the same group of authors?
3. Do authors mainly reference papers with other authors with whom they've co-authored papers (including themselves)?
4. What is the distribution in time for citations of papers in general, and for papers with the highest number of citations? Do they differ?
5. Is there a correlation between citations and the number of keywords that papers share? I.e. papers which share the same subject cite each other more often.
6. For the most-cited papers (define your own cutoff), is the correlation in shared keywords between them and the papers that cite them different from (5) ?
7. Do papers mostly reference other papers that were written in the same language?
8. What is the most cited paper?
9. Who is the most cited author?
10. Which author(s) has/have the highest h-inex?

* * *
## Getting Started
This project is divided into three different sections:
1. Process the data using multiprocessing
2. Create graphs
3. Analyze data using pyspark and networkx

### Processing the data
The PubMed data (XML format) was processed in parallel using the built-in python library [multiprocessing](https://docs.python.org/3/library/multiprocessing.html). The [pubmed_parser](src/pubmed_parser/) and [parallel_computing](src/parallel_computing/) directories contain the code for parsing the data and setting up a server and client for processing the data in parallel. 

To [main.py](src/main.py) script can be used to process the data. It requires the user to specify a few arguments:

<details>
  <summary>main.py arguments</summary>

```bash
  -h, --help          show help message
  -d, --data_dir      Location of the PubMed data (XML)
  -o, --out-dir       The output directory to store the results - 
                      required if server mode is selected
  -p, --port_number   The port number
  --host              Host used
  -n                  Number of CPU used per client
  -s                  Server mode, can not be used together with -c
  -c                  Client mode, can not be used together with -s
```

</details>
> NOTE: to reproduce this help message run: python src/main.py -h
In order to process the data a network([star topology](https://www.computerhope.com/jargon/s/startopo.htm)) needs to be set up where one computer acts as the server that divides the work and one or multiple computers act as client that take work from the server and process the data. 

Therefore, the [main.py](src/main.py) needs to be ran atleast twice(once as a server and once as a client). To specify wheter the script will act as a server or client can be indicated with the `-s` and `-c` flags respectively.

Server:
```bash
python src/main.py -d /data/dataprocessing/NCBI/PubMed/ -p 4235 --host assemblix2019 -s -o /commons/dsls/dsph/2022/test/
```

Client: 
```bash
python src/main.py -p 4235 --host assemblix2019 -c -n 2
```

* * *
## Results

![graph_author](figures/common_authors_graph.png)

* * *
## Installation


* * *
## Requirements

* * *
## Author

Stijn Arends - contact me via [mail](mailto:s.arends@st.hanze.nl).

* * *
## License
This project contains an [MIT license](../LICENSE)

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/stijn-arends/programming3.svg?style=for-the-badge
[contributors-url]: https://github.com/stijn-arends/programming3/graphs/contributors
[license-shield]: https://img.shields.io/github/license/stijn-arends/programming3.svg?style=for-the-badge
[license-url]: https://github.com/stijn-arends/programming3/blob/master/LICENSE.md
