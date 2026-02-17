#!/bin/bash

# Reformatting Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind_underscoreless.vcf into distinct genotype files (one per condition) for various downstream analyses
# This should be run after generating the covariate-regressed files from the R script but before running GenAbel GRAMMAR

# How many CPUs for this job?
#$ -pe parallel 2
# How much RAM per CPU?
#$ -l h_vmem=8G
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o reformat_genotype_files.out
# Use /bin/bash to execute this script
#$ -S /bin/bash
# Run job from current working directory
#$ -cwd
# Send email when the job begins, ends, aborts, or is suspended
#$ -m beas
#$ -M james.tanshengyi@tuebingen.mpg.de

# (1) Initialize the correct conda environment
source ~/.bashrc
conda activate veQTL_GraVe
cd ../data

# (2) Remove underscores from SNP and sample IDs of the provdied file - some downstream software does not like these...
awk '{gsub(/_/, ""); print}' Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind.vcf > Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind_underscoreless.vcf

# (3) Split the vcf file into control and HS with MAF >0.05 and Missingess < 0.5 based on the sample IDs without underscores
vcftools --vcf Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind_underscoreless.vcf --maf 0.05 --max-missing 0.5 --keep Ctrl_samples.txt --recode --out Dmel_Ctrl_MAF5_Miss50_All
vcftools --vcf Dmel_head_hs_ctrl_Miss80_MAF5_LD8_HWE_1975ind_underscoreless.vcf --maf 0.05 --max-missing 0.5 --keep HS_samples.txt --recode --out Dmel_HS_MAF5_Miss50_All

# (4) For each condition, only include SNPs with MAF >0.05 and Missingess < 0.5 common to BOTH conditions. Both files should have the same number of SNPs in the end.
awk '!/^#/ {print $3}' Dmel_Ctrl_MAF5_Miss50_All.recode.vcf | sort >  Dmel_Ctrl_MAF5_Miss50_All.txt
awk '!/^#/ {print $3}' Dmel_HS_MAF5_Miss50_All.recode.vcf | sort >  Dmel_HS_MAF5_Miss50_All.txt
comm -12 Dmel_Ctrl_MAF5_Miss50_All.txt Dmel_HS_MAF5_Miss50_All.txt > common_snps.txt
vcftools --vcf Dmel_Ctrl_MAF5_Miss50_All.recode.vcf --snps common_snps.txt --recode --out Dmel_Ctrl_MAF5_Miss50
vcftools --vcf Dmel_HS_MAF5_Miss50_All.recode.vcf --snps common_snps.txt --recode --out Dmel_HS_MAF5_Miss50

#rm Dmel_Ctrl_MAF5_Miss50_All.recode.vcf Dmel_HS_MAF5_Miss50_All.recode.vcf Dmel_Ctrl_MAF5_Miss50_All.txt Dmel_HS_MAF5_Miss50_All.txt common_snps.txt

# (5) Turn first into .ped and .map files then into .raw files for regressing out two separate GRMs with GenAbel GRAMMAR
plink --vcf Dmel_Ctrl_MAF5_Miss50.recode.vcf --allow-extra-chr --recode tab --out Dmel_Ctrl_MAF5_Miss50_GRAMMAR
plink --vcf Dmel_HS_MAF5_Miss50.recode.vcf --allow-extra-chr --recode tab --out Dmel_HS_MAF5_Miss50_GRAMMAR
Rscript ../code/2a_GRAMMAR_Raw_Genotype_Format.R 

# (6) Turn into .bed files for eQTL mapping and PCA
plink --vcf Dmel_Ctrl_MAF5_Miss50.recode.vcf --allow-extra-chr --make-bed --out Dmel_Ctrl_final_eqtl
plink --vcf Dmel_HS_MAF5_Miss50.recode.vcf --allow-extra-chr --make-bed --out Dmel_HS_final_eqtl
plink --bfile Dmel_Ctrl_final_eqtl --bmerge Dmel_HS_final_eqtl --allow-extra-chr --make-bed --out Dmel_Ctrl_HS_final_full


# (7) Turn into .vcf.gz and .vcf.gz.tbi files for veQTL mapping

# (7a) cis-veQTL mapping 
bgzip -c Dmel_Ctrl_MAF5_Miss50.recode.vcf > Dmel_Ctrl_final_cis_veqtl.vcf.gz 
bgzip -c Dmel_HS_MAF5_Miss50.recode.vcf > Dmel_HS_final_cis_veqtl.vcf.gz 
tabix -p vcf Dmel_Ctrl_final_cis_veqtl.vcf.gz
tabix -p vcf Dmel_HS_final_cis_veqtl.vcf.gz

# (7b) Trans-veQTL mapping - includes renaming SNP ID (from 1 til 383710 - the index of the last SNP) and chromosome (to 2L) steps
grep '^#' Dmel_Ctrl_MAF5_Miss50.recode.vcf > vcf_header.txt  
grep -v '^#' Dmel_Ctrl_MAF5_Miss50.recode.vcf > vcf_data.txt  
awk 'BEGIN {OFS="\t"} {$2=NR; print}' vcf_data.txt > repositioned_vcf_data.txt  
cat vcf_header.txt repositioned_vcf_data.txt > repositioned_Dmel_Ctrl_MAF5_Miss50.vcf
sed '/^[^#]/s/^[^\t]*/2L/' repositioned_Dmel_Ctrl_MAF5_Miss50.vcf > Dmel_Ctrl_final_trans_veqtl.vcf
bgzip -c Dmel_Ctrl_final_trans_veqtl.vcf > Dmel_Ctrl_final_trans_veqtl.vcf.gz 
tabix -p vcf Dmel_Ctrl_final_trans_veqtl.vcf.gz

grep '^#' Dmel_HS_MAF5_Miss50.recode.vcf > vcf_header.txt  
grep -v '^#' Dmel_HS_MAF5_Miss50.recode.vcf > vcf_data.txt  
awk 'BEGIN {OFS="\t"} {$2=NR; print}' vcf_data.txt > repositioned_vcf_data.txt  
cat vcf_header.txt repositioned_vcf_data.txt > repositioned_Dmel_HS_MAF5_Miss50.vcf
sed '/^[^#]/s/^[^\t]*/2L/' repositioned_Dmel_HS_MAF5_Miss50.vcf > Dmel_HS_final_trans_veqtl.vcf
bgzip -c Dmel_HS_final_trans_veqtl.vcf > Dmel_HS_final_trans_veqtl.vcf.gz 
tabix -p vcf Dmel_HS_final_trans_veqtl.vcf.gz

rm vcf_header.txt vcf_data.txt repositioned_vcf_data.txt repositioned_Dmel_Ctrl_MAF5_Miss50.vcf repositioned_Dmel_HS_MAF5_Miss50.vcf

# (8) Create a file containing SNP IDs, chromosome, recoded positions and major/minor alleles from final SNP .vcf files
awk '!/^#/ {print $1, $2, $3, $4, $5}' Dmel_Ctrl_MAF5_Miss50.recode.vcf > SNPs_original_positions.txt
awk '!/^#/ {print $1, $2, $3, $4, $5}' Dmel_Ctrl_final_trans_veqtl.vcf > SNPs_dummy_positions.txt


