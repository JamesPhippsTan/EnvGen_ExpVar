#!/bin/bash

# map_trans_veqtl_HS.sh - Map veQTLs for the HS dataset

# How many CPUs for this job?
#$ -pe parallel 1
# How much RAM per CPU?
#$ -l h_vmem=50G
# Number of parallel jobs to be run as part of an array job
#$ -t 1:293
#$ -o map_trans_veqtl_HS.out
#$ -e map_trans_veqtl_HS.err
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd

# Initialize the correct conda environment
source ~/.bashrc
conda activate veQTL_GraVe

# Map Ctrl trans-veQTL
mkdir -p ../data/results_veqtl_HS_trans
./veqtl-mapper --bed ../data/HS_phenos_trans_veQTL_mapping.bed --vcf ../data/Dmel_HS_final_trans_veqtl.vcf.gz --perm 50,4 --genes 10 --eqtl ../data/eQTL_for_HS_trans_veqtl.txt --job-number $SGE_TASK_ID --out ../data/results_veqtl_HS_trans/trans_batch_${SGE_TASK_ID}