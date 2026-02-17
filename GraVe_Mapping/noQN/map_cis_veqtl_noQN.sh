#!/bin/bash

# map_cis_veqtl_noqn.sh - Map veQTLs for the control dataset

# How many CPUs for this job?
#$ -pe parallel 2
# How much RAM per CPU?
#$ -l h_vmem=32G
# Number of parallel jobs to be run as part of an array job
#$ -t 1:69
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o map_cis_veqtl_noQN.out
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd

# Initialize the correct conda environment
source ~/.bashrc
conda activate veQTL_GraVe

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tmp/global2/jtanshengyi/miniconda3/envs/veQTL_GraVe/lib/
export LD_LIBRARY_PATH

# Map Ctrl cis-veQTL
cd ..
mkdir -p ../data/noQN
./veqtl-mapper-dros-bin --bed ../data/noQN/Ctrl_phenos_cis_mapping_noQN.bed --vcf ../data/Dmel_Ctrl_final_cis_veqtl.vcf.gz --window 10000 --perm 3,4 --genes 127 --eqtl ../data/eQTL_for_Ctrl_cis_veqtl.txt --job-number $SGE_TASK_ID --out ../data/noQN/Ctrl_cis_batch_${SGE_TASK_ID}

# Map HS cis-veQTL
./veqtl-mapper-dros-bin --bed ../data/noQN/HS_phenos_cis_mapping_noQN.bed --vcf ../data/Dmel_HS_final_cis_veqtl.vcf.gz --window 10000 --perm 3,4 --genes 127 --eqtl ../data/eQTL_for_HS_cis_veqtl.txt --job-number $SGE_TASK_ID --out ../data/noQN/HS_cis_batch_${SGE_TASK_ID}
