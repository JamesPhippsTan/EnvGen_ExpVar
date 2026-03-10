#!/bin/bash

# Obtains stats for each SNP within each of the two populations - MAF, genotype counts, and missingness
# SNPs with low stats (esp. missingness) may be used to guide selection of the most interesting veQTL
# e.g., disregard those with low missingess

# How many CPUs for this job?
#$ -pe parallel 4
# How much RAM per CPU?
#$ -l h_vmem=4G
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o mapping_genotypes_stats.out
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
cd ../data

# Get genotype counts using the cis-files as they contain the correct SNP and chromosome positions. 
# Allele freqs and missingness can be calculated in R and appended to the significant results.
plink --bfile Dmel_Ctrl_final_eqtl --allow-extra-chr --freqx --out Dmel_Ctrl_geno_counts
plink --bfile Dmel_HS_final_eqtl --allow-extra-chr --freqx --out Dmel_HS_geno_counts

#  Get per-site genotyping rate missingness
vcftools --vcf Dmel_Ctrl_MAF5_Miss50.recode.vcf --missing-site --out Dmel_Ctrl_SiteMissingness
vcftools --vcf Dmel_HS_MAF5_Miss50.recode.vcf --missing-site --out Dmel_HS_SiteMissingness

# Get per-individual genotyping rate missingness
vcftools --vcf Dmel_Ctrl_MAF5_Miss50.recode.vcf --missing-indv --out Dmel_Ctrl_IndivMissingness
vcftools --vcf Dmel_HS_MAF5_Miss50.recode.vcf --missing-indv --out Dmel_HS_IndivMissingness

