#!/bin/bash

forward=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R1_001_BC24EVACXX.filt.fastq 
reverse=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R2_001_BC24EVACXX.filt.fastq 

out_folder=/commons/dsls/dsph/2022/velvet_out_run

mkdir -p ${out_folder}

# Output destination
output_kmr=${out_folder}/kmr_ # output=${out_folder}/kmr_

# Final output
output=output
mkdir -p ${output}
kmr_size_N50=${output}/kmr_size_N50.csv # ${out_folder}/kmr_size_N50.csv

seq 19 2 31 | parallel -j16 "velveth ${output_kmr}{} {} -longPaired -fastq -separate ${forward} ${reverse} && velvetg ${output_kmr}{} && cat ${output_kmr}{}/contigs.fa | python3 assignment4.py -kmr {} >> ${kmr_size_N50}"

highest_N50=0
highest_kmr=1

while IFS=, read -r kmr n50
do
    if [ "$n50" -gt "$highest_N50" ]; then
    highest_N50=${n50}
    highest_kmr=${kmr}
    fi
done < ${kmr_size_N50}

echo "The highest kmr is ${highest_kmr} with an N50 of: ${highest_N50}"

# Move contigs.fa to the output folder
mv ${output_kmr}${highest_kmr}/contigs.fa ${output}


vals=($(seq 19 2 31))

for kmr in ${vals[@]}
do
    # If you don't want to remove the folder of the best kmr size:
    # ----
    # if [ "$kmr" -ne "$highest_kmr" ]; then
    #     echo "Remove: ${kmr}"
    #     rm -r ${output}${kmr}
    # fi
    rm -r ${output_kmr}${kmr}
done
