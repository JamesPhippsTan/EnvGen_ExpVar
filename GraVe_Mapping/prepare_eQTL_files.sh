#!/bin/bash

# Prepare a file of eQTLs so that their effects will be controlled for during a veQTL mapper run

# How many CPUs for this job?
#$ -pe parallel 1
# How much RAM per CPU?
#$ -l h_vmem=8G
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o prepare_eQTL_file.out
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd
# Send email when the job begins, ends, aborts, or is suspended
#$ -m beas
#$ -M james.tanshengyi@tuebingen.mpg.de

# Initialize the correct conda environment
source ~/.bashrc
conda activate veQTL_GraVe

# (2) Execute R script to create the eQTL file
cd ../code
Rscript 4_eQTL_prep_for_veQTL.R