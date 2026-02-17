#!/bin/bash

# Filter the top trans-veQTL SNPs

# How many CPUs for this job?
#$ -pe parallel 10
# How much RAM per CPU?
#$ -l h_vmem=40G
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o filter_top_trans_veQTL.out
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd
# Send email when the job begins, ends, aborts, or is suspended
#$ -m beas
#$ -M james.tanshengyi@tuebingen.mpg.de

# We need to find the critical highest p-value for the BH correction procedure to get adjusted pvals of below 0.05
# With 3 billion tests this critical value must be much much lower than 0.05
# 0.0005 is a third try to snipe this critical value

# Ctrl
#cd /tmp/global2/jtanshengyi/veQTL_GraVe/data/results_veqtl_Ctrl_trans
#awk 'FNR>1||NR==1' trans_batch* > Ctrl_trans_veQTL
#(head -1 Ctrl_trans_veQTL; tail -n+2 Ctrl_trans_veQTL |  awk '$7 < 0.0005') > Ctrl_trans_veQTL.0_0005.pval.txt
#(head -1 Ctrl_trans_veQTL; tail -n+2 Ctrl_trans_veQTL |  awk '$8 < 0.0005') > Ctrl_trans_veQTL.0_0005.pvalperm.txt

# HS
cd /tmp/global2/jtanshengyi/veQTL_GraVe/data/results_veqtl_HS_trans
awk 'FNR>1||NR==1' trans_batch* > HS_trans_veQTL
(head -1 HS_trans_veQTL; tail -n+2 HS_trans_veQTL |  awk '$7 < 0.0005') > HS_trans_veQTL.0_0005.pval.txt
(head -1 HS_trans_veQTL; tail -n+2 HS_trans_veQTL |  awk '$8 < 0.0005') > HS_trans_veQTL.0_0005.pvalperm.txt


