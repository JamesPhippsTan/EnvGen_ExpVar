#!/bin/bash

# map_eqtl_CPU_Ctrl.sh - Map cis eQTLs Ctrl without quantile normalisation
# Run on a node with GPUs 

# How many CPUs for this job?
# How much RAM per CPU?
#$ -pe parallel 2
#$ -l h_vmem=50G
#$ -l 'hostname=node524|node525|node526|node529'

# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o map_eqtl_CPU_Ctrl_noQN.out
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd

# Initialize the correct conda environment
source ~/.bashrc
conda activate tensorqtl

# Configure the environment to run the cluster's cuda 12.4 installation
. /tmp/global2/maan/hello-cuda-12.4/env

# Map eQTL
cd ../../data/noQN
mkdir -p ./Ctrl_cis_eqtl_result_noQN/
python3 -m tensorqtl ../Dmel_Ctrl_final_eqtl Ctrl_phenos_cis_mapping_noQN.bed Ctrl_noQN --mode cis --window 10000 --seed 4 -o ./Ctrl_cis_eqtl_result_noQN