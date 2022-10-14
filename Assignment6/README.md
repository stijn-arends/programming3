# Assignment 6


## Getting Started

Server:
```bash
python src/main.py -d /data/dataprocessing/NCBI/PubMed/ -p 4235 --host assemblix2019 -s -o /commons/dsls/dsph/2022/test/
```

Client: 
```bash
python src/main.py -p 4235 --host assemblix2019 -c -n 2
```