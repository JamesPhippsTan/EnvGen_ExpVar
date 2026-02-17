#!/bin/bash

# map_eqtl_GPU_Ctrl.sh - Map cis and trans Ctrl eQTLs
# Run on a node with GPUs 

# Run on node with at least 1 GPU with 20GB of memory
#$ -l gpumem=20G 
#$ -l gpu=1
#$ -l h_vmem=50G

# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o map_eqtl_GPU_Ctrl.out
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

mkdir -p ./Ctrl_trans_eqtl_result
python3 -m tensorqtl Dmel_Ctrl_final_eqtl Ctrl_phenos_trans_eQTL_mapping.bed Ctrl --mode trans --output_text --pval_threshold 0.005 --seed 4 -o ./Ctrl_trans_eqtl_result

mkdir -p ./Ctrl_cis_eqtl_result
python3 -m tensorqtl Dmel_Ctrl_final_eqtl Ctrl_phenos_cis_mapping.bed Ctrl --mode cis --window 10000 --seed 4 -o ./Ctrl_cis_eqtl_result

