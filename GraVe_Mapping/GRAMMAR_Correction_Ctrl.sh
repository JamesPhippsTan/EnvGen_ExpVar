#!/bin/bash

# Applying the Grammar Correction R Code on the cluster and generate the phenotype files

# How many CPUs for each parallel job?
#$ -pe parallel 4
# How much RAM per CPU?
#$ -l h_vmem=4G
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o GRAMMAR_Correction_Ctrl.out
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd
# Send email when the job begins, ends, aborts, or is suspended
#$ -m beas
#$ -M james.tanshengyi@tuebingen.mpg.de
# Number of parallel jobs to be run as part of an array job - 40 parallel jobs
#$ -t 1:40

# Initialize the correct conda environment
source ~/.bashrc
conda activate veQTL_GraVe

# Run the Rscript to get GRAMMAR-corrected phenotypes
Rscript 2b_Phenotypes_GRAMMAR_Corrected_Ctrl.R ${SGE_TASK_ID}