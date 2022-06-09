#!/bin/bash

forward=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R1_001_BC24EVACXX.filt.fastq 
reverse=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R2_001_BC24EVACXX.filt.fastq 

out_folder=/commons/dsls/dsph/2022/velvet_out_run_full/output

mkdir -p ${out_folder}

# Output destination
output=${out_folder}/kmr_

kmr_size_N50=${out_folder}/kmr_size_N50.csv


# Works
# -------

# seq 27 2 31 | parallel -j16 velveth $output{} {} -longPaired -fastq -separate $forward $reverse

# seq 27 2 31 | parallel velvetg $output{}

# Best option
seq 5 2 31 | parallel -j16 "velveth ${output}{} {} -longPaired -fastq -separate ${forward} ${reverse} && velvetg ${output}{} && cat ${output}{}/contigs.fa | python3 assignment4.py -kmr {} >> ${kmr_size_N50}"


# Hardcoded paths
# seq 27 2 31 | parallel -j16 'velveth /commons/dsls/dsph/2022/test2/Stijn_{} {} -longPaired -fastq -separate /data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R1_001_BC24EVACXX.filt.fastq /data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R2_001_BC24EVACXX.filt.fastq && velvetg /commons/dsls/dsph/2022/test2/Stijn_{}'

# -------

# Find the best kmr size
# ----------------
highest_N50=0
highest_kmr=1

while IFS=, read -r kmr n50
do
    # echo "I got:$kmr|$n50"
    if [ "$n50" -gt "$highest_N50" ]; then
    # echo "${n50} > ${highest_N50}"
    highest_N50=${n50}
    highest_kmr=${kmr}
    fi
done < ${kmr_size_N50}

echo "The highest kmr is ${highest_kmr} with an N50 of: ${highest_N50}"

# Move contigs.fa to the output folder
mv ${output}${highest_kmr}/contigs.fa ${out_folder}


vals=($(seq 5 2 31))

for kmr in ${vals[@]}
do
    # If you don't want to remove the folder of the best kmr size:
    # ----
    # if [ "$kmr" -ne "$highest_kmr" ]; then
    #     echo "Remove: ${kmr}"
    #     rm -r ${output}${kmr}
    # fi
    rm -r ${output}${kmr}
done
