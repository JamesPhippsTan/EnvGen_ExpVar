#!/bin/bash

# Run SnpEff on my SNPs

# How many CPUs for this job?
#$ -pe parallel 2
# How much RAM per CPU?
#$ -l h_vmem=8G
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o SnpEff.out
# Use /bin/bash to execute this script
#$ -S /bin/bash
# Run job from current working directory
#$ -cwd
# Send email when the job begins, ends, aborts, or is suspended
#$ -m beas
#$ -M james.tanshengyi@tuebingen.mpg.de

# (1) Initialize the correct conda environment
source ~/.bashrc

# (2) Get SNPeff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip

# (3) Get a compatible java (openjdk) version (21.0.6 or above) and create a conda environment called java for this purpose
conda search -f openjdk # search available versions
conda create -n java openjdk=21.0.6=h68779a4_0
conda activate java 

# (5) Eliminate the sample genotypes and change the chromosome '23' annotation to 'X'
cd ../data
awk '{if($0 ~ /^#/){print $0} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind_underscoreless.vcf > SNP_Set.vcf
awk 'BEGIN {FS=OFS="\t"} {if ($1 == "23") $1 = "X"; print}' SNP_Set.vcf > SNP_Set_X.vcf


# (6) Check the drosophila databases available to find the most appropriate genome SNP annotation
java -jar ../code/snpEff/snpEff.jar databases | grep 'melanogaster'
# BDGP6.28.99 works

# (5) Get the SNP annotations
java -Xmx8g -jar ../code/snpEff/snpEff.jar -v -csvStats SnpEff_Annotated_SNP_Set.csv -stats SnpEff_Annotated_SNP_Set.html BDGP6.28.99 SNP_Set_X.vcf > SnpEff_Annotated_SNP_Set.vcf
awk '!/^#/ {print $1, $2, $3, $4, $5, $6, $7, $8}' SnpEff_Annotated_SNP_Set.vcf > SnpEff_Annotated_SNP_Set.txt
rm SnpEff_Annotated_SNP_Set.vcf SNP_Set_X.vcf



