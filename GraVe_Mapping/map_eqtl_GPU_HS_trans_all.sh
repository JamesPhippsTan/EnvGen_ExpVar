#!/bin/bash

# map_eqtl_GPU_HS_trans_all.sh - Map all trans Ctrl eQTLs
# Run on a node with GPUs 

# Run on node with at least 1 GPU with 20GB of memory
#$ -l gpumem=20G 
#$ -l gpu=1
#$ -l h_vmem=50G

# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o map_eqtl_GPU_HS_trans_all.out
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd
# Send email when the job begins, ends, aborts, or is suspended
#$ -m eas
#$ -M james.tanshengyi@tuebingen.mpg.de

# Initialize the correct conda environment
source ~/.bashrc
conda activate tensorqtl

# Configure the environment to run the cluster's cuda 12.4 installation
. /tmp/global2/maan/hello-cuda-12.4/env

# Map eQTL
cd ../data

mkdir -p ./HS_trans_eqtl_result_all
python3 -m tensorqtl Dmel_HS_final_eqtl HS_phenos_trans_eQTL_mapping.bed Ctrl --mode trans --output_text --pval_threshold 1 --seed 4 -o ./HS_trans_eqtl_result_all
