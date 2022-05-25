# Assignment 3 of the porgramming 3 course
* * *

* Date: 21-5-2022
* Author: Stijn Arends
* Assignment: [link](https://bioinf.nl/~martijn/master/programming3/assignment3.html)

* * *
## Assignment

You need to submit a bash script called assignment3.sh in the Assignment3 folder of your programming3 repository on GitHub. 
This script tests the blastp program with different "-numthreads" on the query file "MCRA.faa" (an Amino Acid file of the McrA protein involved in biogas production). 
This script makes an output directory called output when run, and places the run times for each blastp run in a file called timings.txt in the output directory. 
The timings should be based on the Linux "time" command (/usr/bin/time; make sure you're not using the built-in shell version, 
use the whole command name with path and all).

Your assignment is to run the program with -num_threads values from 1 to 16 and capture the time spent running in seconds

```bash
export BLASTDB=/local-fs/datasets/
blastp -query MCRA.faa -db refseq_protein/refseq_protein -num_threads 1 -outfmt 6 >> blastoutput.txt
```

The output should be a file called timings.txt and a PNG format matplotlib graph showing number-of-threads on the x axis and time-taken (s) on the y axis, called timings.png. (Note: so you need to run a small python script as well to generate the graph).

## Running the script

```bash
bash ./assignment3.sh
```

* * *
## Packages

| Package                                                           | Version        |
| ----------------------------------------------------------------- | :------------: |
| [matplotlib](https://matplotlib.org/)                             | `3.3.4`        |
| [pandas](https://pandas.pydata.org/)                              | `1.1.5`        |
