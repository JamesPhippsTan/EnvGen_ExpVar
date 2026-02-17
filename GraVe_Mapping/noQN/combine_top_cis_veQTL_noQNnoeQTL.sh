#!/bin/bash

# Combine the top cis-veQTL SNPs (no eQTL tag)

# How many CPUs for this job?
#$ -pe parallel 2
# How much RAM per CPU?
#$ -l h_vmem=10G
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o Combine_top_cis_veQTL_noQN.out
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd
# Send email when the job begins, ends, aborts, or is suspended
#$ -m beas
#$ -M james.tanshengyi@tuebingen.mpg.de

# Ctrl
cd /tmp/global2/jtanshengyi/veQTL_GraVe/data/noQN
awk 'FNR>1||NR==1' Ctrl_cis_batch_noeQTL* > Ctrl_cis_veQTL_noQNnoeQTL.txt
awk 'FNR>1||NR==1' HS_cis_batch_noeQTL* > HS_cis_veQTL_noQNnoeQTL.txt

