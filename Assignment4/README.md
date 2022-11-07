# Assignment 4 of the porgramming 3 course
* * *

* Date: 01-6-2022
* Author: Stijn Arends
* Assignment: [link](https://bioinf.nl/~martijn/master/programming3/assignment4.html)

* * *
## Assignment

The deliverable is a Bash script called assignment4.sh and a Python script called assignment4.py in the Assignment4/ directory of your programming3 GitHub repository.

The NGS data you should assemble is located at `/data/dataprocessing/MinIONData/MG5267/` on the assemblix servers.

Your bash script should control the job; it should run `velveth` and `velvetg`, using GNU parallel of course to parallelize trying different kmer sizes. Your python script should examine the produced assemblies (they will be in FASTA format) and calculate the N50. Your python script should read any number of contigs from `sys.stdin` and write the N50 of those contigs to `sys.stdout` (as a simple number). GNU Parallel controls how many lines are sent into your python script using its chunking facilities, and using up to 16 jobs in parallel.

Your final output should be located in an `output/` directory and should be the "best" assembly FASTA file as measured by N50, as well as a CSV file containing the kmer parameters tried and the N50 of the resulting assembly. (Be sure to delete all the other assembly files produced by your script except the best one)


## Running the script

```bash
bash ./assignment4.sh
```

* * *
## Packages

| Package                                                           | Version        |
| ----------------------------------------------------------------- | :------------: |
| [numpy](https://numpy.org/)                                       | `1.23.3`       |

