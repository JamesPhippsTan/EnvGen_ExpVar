#!/bin/bash

# map_trans_veqtl_Ctrl.sh - Map veQTLs for the control dataset

# How many CPUs for this job?
#$ -pe parallel 1
# How much RAM per CPU?
#$ -l h_vmem=48G
# Number of parallel jobs to be run as part of an array job
#$ -t 1:438
#$ -o map_trans_veqtl_Ctrl.out
#$ -e map_trans_veqtl_Ctrl.err
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd
#$ -m a
#$ -M james.tanshengyi@tuebingen.mpg.de

# Initialize the correct conda environment
source ~/.bashrc
conda activate veQTL_GraVe

# Map Ctrl trans-veQTL
mkdir -p ../data/results_veqtl_Ctrl_trans
./veqtl-mapper --bed ../data/Ctrl_phenos_trans_veQTL_mapping.bed --vcf ../data/Dmel_Ctrl_final_trans_veqtl.vcf.gz --perm 2,4 --genes 10 --eqtl ../data/eQTL_for_Ctrl_trans_veqtl.txt --job-number $SGE_TASK_ID --out ../data/results_veqtl_Ctrl_trans/trans_batch_${SGE_TASK_ID}