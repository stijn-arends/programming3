#!/bin/bash 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stijnarends@live.nl
#SBATCH --time 2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=StijnAssg3
#SBATCH --partition=assemblix 

export cpus=16
export BLASTDB=/local-fs/datasets/

# echo "Number of cpus:${cpus}"

#for i in {1..16}; do echo "Cpu: ${i}"; done

mkdir output

# for i in $( seq 1 $cpus )
# do 
#     /usr/bin/time -a --output output/timings.txt -f %e blastp -query MRCA.faa -db ${BLASTDB}refseq_protein/refseq_protein -num_threads $i -outfmt 6 >> blastoutput.txt
#     # awk 'BEGIN{ printf " $i" >> "output/timings.txt" }'
#     sed -i "$ s/^/$i /" output/timings.txt
# done

for i in {1..16}
do 
    # echo "Processing using ${i} cores"
    /usr/bin/time -a --output output/timings_check.txt -f %e blastp -query MRCA.faa -db ${BLASTDB}refseq_protein/refseq_protein -num_threads $i -outfmt 6 >> blastoutput.txt 
    sed -i "$ s/^/$i /" output/timings_check.txt
done

python3 plot_time_benefit.py output/timings_check.txt