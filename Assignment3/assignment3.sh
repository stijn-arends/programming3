#!/bin/bash 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stijnarends@live.nl
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=BlastpStijn
#SBATCH --partition=assemblix

# CPUs
cpus=16

# Get the relative directory of the script
script_relative_dir=$(dirname "${BASH_SOURCE[0]}") 
echo $script_relative_dir

# Output 
out_folder=${script_relative_dir}/output
out_time=${out_folder}/timings.txt
out_blast=blastoutput.txt

# Input
query_file=MCRA.faa
export BLASTDB=/local-fs/datasets/

echo "Output folder: ${out_folder}"

if [ ! -f "${query_file}" ]; then
    echo "Invalid query file. File ${query_file} does not exist. Make sure that you are in the correct directory."
    exit
fi

if [ ! -d "${out_folder}" ]; then
    mkdir ${out_folder}
    echo "Output does not exists yet, creating..."
fi


for i in $( seq 16 $cpus )
do
    echo "Iteration: $i"
    /usr/bin/time -a --output ${out_time} -f "$i %e" blastp -query ${query_file} -db refseq_protein/refseq_protein -num_threads $i -outfmt 6 >> ${out_blast}
done

python3 plot_time_benefit.py ${out_time}