# Investigating eQTL and veQTL 
# Part 2: Null distribution generation for derived allele effects
# This script creates the selection criteria files for non-eQTL
# Huiting's bash pipeline runs the subsampling to create null distributions
# The resultant summary statistics are then plugged back into this file
# to generate plots

# Last Updated: 22/9/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls(all = T))

library(reshape)
library(tidyr)
library(dplyr)
library(stats)
library(EnvStats)
library(data.table)
library(ggplot2)
library(UpSetR)
library(tibble)
library(ggvenn)
library(gridExtra)
library(patchwork)
library(purrr)

# Load saved environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
load(file='6_Investigating_eQTL_and_veQTL_Allele_Age2.RData')

# Functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('gProfiler_Functions.R')
source('Quantile_Functions.R')
source('Variability_Functions.R')
source('QTL_Analysis_Functions.R') 

#######################
##### Datasets ########
#######################

# Continues with objects saved in the following .RData files
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
load(file='5_Investigating_eQTL_and_veQTL.RData')
load(file='6_Investigating_eQTL_and_veQTL_Allele_Age1.RData')

# MAFs for each SNP in each condition
View(MAF)
# Narrow them down to the MAFs that are present in the file
MAF <- MAF[SNPs_position_map$SNP,]

# Allele ages from previous dataframe

# Get all cis-SNPs with age
all_cis_SNPs <- unique(Ctrl_cis_veQTL$SNP)
length(all_cis_SNPs)
length(intersect(all_cis_SNPs,unique(SNP_allele_age$SNP)))
# 104637

# Get all (trans) SNPs with age
all_trans_SNPs <- SNPs_position_map$SNP
length(all_trans_SNPs)
length(intersect(all_trans_SNPs,unique(SNP_allele_age$SNP)))
# 121683


##################################################################
# (3) ....more than any equivalent set of non-eQTL or non-veQTL? #
##################################################################

# Used to create non-veQTL populations of equal size and MAFs of the current sets
#(Note: Not using # of genes this time round for time concerns)

#########################
# (3a) ....non-cis-eQTL #
#########################

# Save criteria to the following directory
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Derived_Allele_Increased_Fraction_Null")

# MAF Criteria.txt
Ctrl_FDI_SNP_MAF <- MAF[Ctrl_cis_per_eQTL_FDI$SNP,'Ctrl_MAF']
Ctrl_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  Ctrl_MAF_Criteria[quintile,'Min']=quantile(Ctrl_FDI_SNP_MAF,min)
  Ctrl_MAF_Criteria[quintile,'Max']=quantile(Ctrl_FDI_SNP_MAF,max)
  Ctrl_MAF_Criteria[quintile,'NumSNPstoSample']=round(length(Ctrl_FDI_SNP_MAF)/5)
}
View(Ctrl_MAF_Criteria)
write.table(Ctrl_MAF_Criteria,'Ctrl_cis-eQTL_MAF_Criteria.txt',row.names = F,col.names = T,quote = F)
# Repeat for HS
HS_FDI_SNP_MAF <- MAF[HS_cis_per_eQTL_FDI$SNP,'HS_MAF']
HS_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  HS_MAF_Criteria[quintile,'Min']=quantile(HS_FDI_SNP_MAF,min)
  HS_MAF_Criteria[quintile,'Max']=quantile(HS_FDI_SNP_MAF,max)
  HS_MAF_Criteria[quintile,'NumSNPstoSample']=round(length(HS_FDI_SNP_MAF)/5)
}
View(HS_MAF_Criteria)
write.table(HS_MAF_Criteria,'HS_cis-eQTL_MAF_Criteria.txt',row.names = F,col.names = T,quote = F)

# non-eQTL_SNPs.txt
# Needs to have age and NOT be a eQTL
non_Ctrl_cis_eQTL<- intersect(setdiff(all_cis_SNPs,Ctrl_cis_per_eQTL_FDI$SNP),unique(SNP_allele_age$SNP))
length(non_Ctrl_cis_eQTL)
# 97K of these 
# Obtain the MAFs
non_Ctrl_cis_eQTL_by_MAF <- data.frame(SNP=non_Ctrl_cis_eQTL,MAF=MAF[non_Ctrl_cis_eQTL,'Ctrl_MAF'])
View(non_Ctrl_cis_eQTL_by_MAF)
write.table(non_Ctrl_cis_eQTL_by_MAF,'Ctrl_non-cis-eQTL_SNPs.txt',row.names = F,col.names = T,quote = F)
# Repeat for HS
non_HS_cis_eQTL<- intersect(setdiff(all_cis_SNPs,HS_cis_per_eQTL_FDI$SNP),unique(SNP_allele_age$SNP))
length(non_HS_cis_eQTL)
# 100K of these 
# Obtain the MAFs
non_HS_cis_eQTL_by_MAF <- data.frame(SNP=non_HS_cis_eQTL,MAF=MAF[non_HS_cis_eQTL,'HS_MAF'])
View(non_HS_cis_eQTL_by_MAF)
write.table(non_HS_cis_eQTL_by_MAF,'HS_non-cis-eQTL_SNPs.txt',row.names = F,col.names = T,quote = F)

###########################
# (3b) ....non-trans-eQTL #
###########################

# The FDD and FDI result was quite striking
# dataframe_name	Number_of_SNPs_w_Age	FDI_all	FDD_all	In_between
# Ctrl_trans_eqtl	92604	25.68	29.37778066	44.94406289
# HS_trans_eqtl	66997	32.31	36.39565951	31.29095333
# Most of the SNPs that have allele age (121K total) are eQTL (90K and 60K)

# The fraction is still large when considering only those with an FDI of 1 or 0
# All_FDD_or_FDI	in_between	remaining_SNPs_w/_age
# 50985.7072	41620	27394.2928
# 46030.7307	20964	53005.2693

# One idea to get a statistic is to run downsampling FDI FDD 
# On both the eQTL and the non-eQTL
# The All_FDD_or_FDI SNPs - should recapitulate the FDD and FDI - 25:29 ~= 50:50
# The non-eQTL (remaining SNPs w/ age) - should be close to 50:50
# Given the above numbers, 1000 SNPs per sample seems reasonable
# Prepare the big trans-eQTL files from which to sample from based on the above sampling space

# Save criteria to the following directory
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Derived_Allele_Increased_Fraction_Null")

# Let us define the eQTL and non-eQTL SNPs first because this will make this easier in the next steps
# Include those that are intermediate as well
Ctrl_trans_eQTL_FDI <- Ctrl_trans_per_eQTL_FDI$SNP
length(Ctrl_trans_eQTL_FDI)
HS_trans_eQTL_FDI <- HS_trans_per_eQTL_FDI$SNP
length(HS_trans_eQTL_FDI)
# No to be diluted by those that are intermediate...
Ctrl_trans_non_eQTL_aged <- intersect(setdiff(all_trans_SNPs,Ctrl_trans_per_eQTL_FDI$SNP),SNP_allele_age$SNP)
length(Ctrl_trans_non_eQTL_aged)
HS_trans_non_eQTL_aged <- intersect(setdiff(all_trans_SNPs,HS_trans_per_eQTL_FDI$SNP),SNP_allele_age$SNP)
length(HS_trans_non_eQTL_aged)

# (1) MAF Criteria - based off of the eQTL SNPs
# Ctrl trans-eQTL MAF Criteria.txt
Ctrl_FDI_SNP_MAF <- MAF[Ctrl_trans_eQTL_FDI,'Ctrl_MAF']
Ctrl_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  Ctrl_MAF_Criteria[quintile,'Min']=quantile(Ctrl_FDI_SNP_MAF,min)
  Ctrl_MAF_Criteria[quintile,'Max']=quantile(Ctrl_FDI_SNP_MAF,max)
  Ctrl_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
View(Ctrl_MAF_Criteria)
write.table(Ctrl_MAF_Criteria,'Ctrl_trans-eQTL_MAF_Criteria_v2.txt',row.names = F,col.names = T,quote = F)

# Ctrl trans-eQTL MAF Criteria.txt
HS_FDI_SNP_MAF <- MAF[HS_trans_eQTL_FDI,'HS_MAF']
HS_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  HS_MAF_Criteria[quintile,'Min']=quantile(HS_FDI_SNP_MAF,min)
  HS_MAF_Criteria[quintile,'Max']=quantile(HS_FDI_SNP_MAF,max)
  HS_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
View(HS_MAF_Criteria)
write.table(HS_MAF_Criteria,'HS_trans-eQTL_MAF_Criteria_v2.txt',row.names = F,col.names = T,quote = F)

# (2) Defining subsampling sets
# eQTL_SNPs and their MAFs
Ctrl_trans_eQTL_by_MAF <- data.frame(SNP=Ctrl_trans_eQTL_FDI,MAF=MAF[Ctrl_trans_eQTL_FDI,'Ctrl_MAF'])
View(Ctrl_trans_eQTL_by_MAF)
write.table(Ctrl_trans_eQTL_by_MAF,'Ctrl_trans-eQTL_downsubsample_SNPs_v2.txt',row.names = F,col.names = T,quote = F)
# Repeat for HS
HS_trans_eQTL_by_MAF <- data.frame(SNP=HS_trans_eQTL_FDI,MAF=MAF[HS_trans_eQTL_FDI,'HS_MAF'])
View(HS_trans_eQTL_by_MAF)
write.table(HS_trans_eQTL_by_MAF,'HS_trans-eQTL_downsubsample_SNPs_v2.txt',row.names = F,col.names = T,quote = F)

# non-eQTL_SNPs and their MAFs
non_Ctrl_trans_eQTL_by_MAF <- data.frame(SNP=Ctrl_trans_non_eQTL_aged,MAF=MAF[Ctrl_trans_non_eQTL_aged,'Ctrl_MAF'])
View(non_Ctrl_trans_eQTL_by_MAF)
write.table(non_Ctrl_trans_eQTL_by_MAF,'Ctrl_non-trans-eQTL_downsubsample_SNPs_v2.txt',row.names = F,col.names = T,quote = F)
# Repeat for HS
non_HS_trans_eQTL_by_MAF <- data.frame(SNP=HS_trans_non_eQTL_aged,MAF=MAF[HS_trans_non_eQTL_aged,'HS_MAF'])
View(non_HS_trans_eQTL_by_MAF)
write.table(non_HS_trans_eQTL_by_MAF,'HS_non-trans-eQTL_downsubsample_SNPs_v2.txt',row.names = F,col.names = T,quote = F)

# (3) Defining datasets to downsample from
# The significant trans-eQTL files have been created previously. I do not need to create it again. I just need to push these to Huiting.

# The non-significant trans-eQTL files have not been created. 
# Creating the entire set is not feasible with the way tensorQTL is coded and my cluster testing
# One way is to narrow down subset the vcf files to just those
# SNPs within 'HS_trans_non_eQTL_aged' and 'Ctrl_trans_non_eQTL_aged'
# Assuming the full file for 400K SNPs is 200GB, that means for the 50K SNPs is 25GB
# This should be more manageable...according to my theory. Let us try this out in practice...
write.table(Ctrl_trans_non_eQTL_aged,'Ctrl_non-trans-eQTL_aged_SNPs.txt',row.names = F,col.names = F,quote = F)
write.table(HS_trans_non_eQTL_aged,'HS_non-trans-eQTL_aged_SNPs.txt',row.names = F,col.names = F,quote = F)

##########################
# (3c) ....non-cis-veQTL #
##########################

# 4 SNPs and 19 SNPs 
# SNP set enrichments are meaningless


############################
# (3d) ....non-trans-veQTL #
############################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Derived_Allele_Increased_Fraction_Null\\trans_veQTL_v2/")

# Let us define the veQTL and non-veQTL SNPs first because this will make this easier in the next steps
Ctrl_trans_veQTL_FDI <- Ctrl_trans_per_veQTL_FDI$SNP
length(Ctrl_trans_veQTL_FDI)
HS_trans_veQTL_FDI <- HS_trans_per_veQTL_FDI$SNP
length(HS_trans_veQTL_FDI)
Ctrl_trans_non_veQTL_aged <- intersect(setdiff(all_trans_SNPs,Ctrl_trans_per_veQTL_FDI$SNP),SNP_allele_age$SNP)
length(Ctrl_trans_non_veQTL_aged)
HS_trans_non_veQTL_aged <- intersect(setdiff(all_trans_SNPs,HS_trans_per_veQTL_FDI$SNP),SNP_allele_age$SNP)
length(HS_trans_non_veQTL_aged)

# (1) MAF Criteria - based off of the veQTL SNPs
# Ctrl trans-veQTL MAF Criteria.txt
Ctrl_FDI_SNP_MAF <- MAF[Ctrl_trans_veQTL_FDI,'Ctrl_MAF']
Ctrl_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  Ctrl_MAF_Criteria[quintile,'Min']=quantile(Ctrl_FDI_SNP_MAF,min)
  Ctrl_MAF_Criteria[quintile,'Max']=quantile(Ctrl_FDI_SNP_MAF,max)
  Ctrl_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
View(Ctrl_MAF_Criteria)
write.table(Ctrl_MAF_Criteria,'Ctrl_trans-veQTL_MAF_Criteria_v2.txt',row.names = F,col.names = T,quote = F)
# Ctrl trans-veQTL MAF Criteria.txt

HS_FDI_SNP_MAF <- MAF[HS_trans_veQTL_FDI,'HS_MAF']
HS_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  HS_MAF_Criteria[quintile,'Min']=quantile(HS_FDI_SNP_MAF,min)
  HS_MAF_Criteria[quintile,'Max']=quantile(HS_FDI_SNP_MAF,max)
  HS_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
View(HS_MAF_Criteria)
write.table(HS_MAF_Criteria,'HS_trans-veQTL_MAF_Criteria_v2.txt',row.names = F,col.names = T,quote = F)


# (2) Defining subsampling sets
# veQTL_SNPs and their MAFs
Ctrl_trans_veQTL_by_MAF <- data.frame(SNP=Ctrl_trans_veQTL_FDI,MAF=MAF[Ctrl_trans_veQTL_FDI,'Ctrl_MAF'])
View(Ctrl_trans_veQTL_by_MAF)
write.table(Ctrl_trans_veQTL_by_MAF,'Ctrl_trans-veQTL_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F)
# Repeat for HS

HS_trans_veQTL_by_MAF <- data.frame(SNP=HS_trans_veQTL_FDI,MAF=MAF[HS_trans_veQTL_FDI,'HS_MAF'])
View(HS_trans_veQTL_by_MAF)
write.table(HS_trans_veQTL_by_MAF,'HS_trans-veQTL_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F)

# non-veQTL_SNPs and their MAFs
non_Ctrl_trans_veQTL_by_MAF <- data.frame(SNP=Ctrl_trans_non_veQTL_aged,MAF=MAF[Ctrl_trans_non_veQTL_aged,'Ctrl_MAF'])
View(non_Ctrl_trans_veQTL_by_MAF)
write.table(non_Ctrl_trans_veQTL_by_MAF,'Ctrl_non-trans-veQTL_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F)
# Repeat for HS
non_HS_trans_veQTL_by_MAF <- data.frame(SNP=HS_trans_non_veQTL_aged,MAF=MAF[HS_trans_non_veQTL_aged,'HS_MAF'])
View(non_HS_trans_veQTL_by_MAF)
write.table(non_HS_trans_veQTL_by_MAF,'HS_non-trans-veQTL_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F)

#################################################
# (4) Looking at the derived allele frequencies #
#################################################

# Convert minor allele frequency into derived allele frequency
MAF$SNP <- rownames(MAF)
MAF_DAF <- merge(MAF,SNP_allele_age[,c('SNP','ANC')],by='SNP')
# Step 1: check if the minor allele is the ancestral allele 
# Step 2: if true, carry out DAF=1-MAF; if false, DAF=MAF
MAF_DAF$Ctrl_DAF <- ifelse(MAF_DAF$Ctrl_minor_allele==MAF_DAF$ANC,1-MAF_DAF$Ctrl_MAF,MAF_DAF$Ctrl_MAF) 
MAF_DAF$HS_DAF <- ifelse(MAF_DAF$HS_minor_allele==MAF_DAF$ANC,1-MAF_DAF$HS_MAF,MAF_DAF$HS_MAF) 
# Check out the MAF and DAF spectrum
hist(MAF_DAF$Ctrl_MAF)
hist(MAF_DAF$HS_MAF)
# Still skewed towards lower-frequency minor alleles
hist(MAF_DAF$Ctrl_DAF)
hist(MAF_DAF$HS_DAF)

# Contains mix of minor alleles AND major alleles
# Still skewed towards lower-frequency derived alleles
MAF_DAF <- column_to_rownames(MAF_DAF,var = 'SNP')
# How many times where the derived allele is the minor allele
table(MAF_DAF$Ctrl_MAF==MAF_DAF$Ctrl_DAF)
table(MAF_DAF$HS_MAF==MAF_DAF$HS_DAF)
# 40K derived is the major allele, 80K derived is the minor allele
# How do I proceed?

# Split by whether the derived allele increases or decreases variability
HS_trans_veQTL_der_inc <- subset(HS_trans_per_veQTL_FDI,FDI==1)
HS_trans_veQTL_der_dec <- subset(HS_trans_per_veQTL_FDI,FDI==0)
HS_trans_veQTL_der_both <- subset(HS_trans_per_veQTL_FDI,FDI<1 & FDI>0)

# Check if true for trans-eQTL
HS_trans_eQTL_der_inc <- subset(HS_trans_per_eQTL_FDI,FDI==1)
HS_trans_eQTL_der_dec <- subset(HS_trans_per_eQTL_FDI,FDI==0)
HS_trans_eQTL_der_both <- subset(HS_trans_per_eQTL_FDI,FDI<1 & FDI>0)

  # Repeat for control
Ctrl_trans_veQTL_der_inc <- subset(Ctrl_trans_per_veQTL_FDI,FDI==1)
Ctrl_trans_veQTL_der_dec <- subset(Ctrl_trans_per_veQTL_FDI,FDI==0)
Ctrl_trans_veQTL_der_both <- subset(Ctrl_trans_per_veQTL_FDI,FDI<1 & FDI>0)

Ctrl_trans_eQTL_der_inc <- subset(Ctrl_trans_per_eQTL_FDI,FDI==1)
Ctrl_trans_eQTL_der_dec <- subset(Ctrl_trans_per_eQTL_FDI,FDI==0)
Ctrl_trans_eQTL_der_both <- subset(Ctrl_trans_per_eQTL_FDI,FDI<1 & FDI>0)

# Function to make one histogram per category
# Colours based on previous code 
make_daf_histogram <- function(df, df_name, MAF_DAF) {
  # Define fixed colours
  my_colors <- c(
    eQTL_der_inc = "#2B5D29",  
    eQTL_der_dec = "#77BB75",  
    veQTL_der_inc = "#611BB8", 
    veQTL_der_dec = "#CFB2F3"  
  )
  
  # Parse condition (Ctrl/HS) and category (eQTL_der_inc etc.)
  condition <- sub("_.*", "", df_name)
  category  <- sub("^(Ctrl|HS)_", "", df_name)
  category  <- sub("^trans_", "", category)   # <-- strip 'trans_'
  daf_col   <- paste0(condition, "_DAF")
  
  # Extract SNPs + DAF values
  snps <- df$SNP
  daf_values <- MAF_DAF[snps, daf_col]
  
  plot_data <- data.frame(DAF = daf_values)
  
  ggplot(plot_data, aes(x = DAF)) +
    geom_histogram(
      bins = 10,
      fill = my_colors[[category]],  # auto-pick colour based on df_name
      color = "black"
    ) +
    labs(
      x = "DAF",
      y = "Count"
    ) +
    theme_classic()
}

# ------------------------------------------------
# Example: collect SNP dataframes into a list
# (replace with your real dataframes)
snp_dfs <- list(
  Ctrl_trans_eQTL_der_inc = Ctrl_trans_eQTL_der_inc,
  Ctrl_trans_eQTL_der_dec = Ctrl_trans_eQTL_der_dec,
  Ctrl_trans_veQTL_der_inc = Ctrl_trans_veQTL_der_inc,
  Ctrl_trans_veQTL_der_dec = Ctrl_trans_veQTL_der_dec,
  HS_trans_eQTL_der_inc = HS_trans_eQTL_der_inc,
  HS_trans_eQTL_der_dec = HS_trans_eQTL_der_dec,
  HS_trans_veQTL_der_inc = HS_trans_veQTL_der_inc,
  HS_trans_veQTL_der_dec = HS_trans_veQTL_der_dec
)

# ------------------------------------------------
# Generate all plots in one step
plots <- imap(snp_dfs, ~ make_daf_histogram(.x, .y, MAF_DAF))

# Show one plot
plots$HS_trans_veQTL_der_inc
plots$HS_trans_veQTL_der_dec

# Save the plots
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\DAF_downsampling")
imap(plots, ~ ggsave(
  filename = paste0(.y, "_DAF_histogram.svg"),
  plot = .x, width = 2, height = 1.5, dpi=300
))

# There is a large skew in towards variability-increasing alleles being high DAF

perform_Wtest <- function(group1_name, group2_name, condition, list_data, maf_df) {
  # Extract the SNPs for each group
  snps1 <- list_data[[group1_name]]$SNP
  snps2 <- list_data[[group2_name]]$SNP
  
  # Select which DAF column to use
  daf_col <- paste0(condition, "_DAF")
  
  # Extract DAF values by merging with MAF_DAF
  daf1 <- maf_df[snps1, daf_col]
  daf2 <- maf_df[snps2, daf_col]
  
  # Perform Wilcoxon ranked sum test
  test_result <- wilcox.test(daf1, daf2)
  
  data.frame(
    comparison = paste(group1_name, "vs", group2_name),
    Wilcoxon_p_value = test_result$p.value,
    median_group1 = median(daf1, na.rm = TRUE),
    n_group1 = length(daf1),
    median_group2 = median(daf2, na.rm = TRUE),
    n_group2 = length(daf2),
    median_diff = median(daf1, na.rm = TRUE) - median(daf2, na.rm = TRUE)
  )
}

# Perform the 4 tests
results <- do.call(rbind, list(
  perform_Wtest("HS_trans_eQTL_der_inc",   "HS_trans_eQTL_der_dec",   "HS",   snp_dfs, MAF_DAF),
  perform_Wtest("Ctrl_trans_eQTL_der_inc", "Ctrl_trans_eQTL_der_dec", "Ctrl", snp_dfs, MAF_DAF),
  perform_Wtest("HS_trans_veQTL_der_inc",  "HS_trans_veQTL_der_dec",  "HS",   snp_dfs, MAF_DAF),
  perform_Wtest("Ctrl_trans_veQTL_der_inc","Ctrl_trans_veQTL_der_dec","Ctrl", snp_dfs, MAF_DAF)
))


setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Derived_Allele_Increased_Fraction/")
write.csv(results,'DAF_derinderdec_Wilcox_tests.csv',row.names = F,col.names = T,quote = F)

##########################################################
# (5) Does minor allele increase or decrease more often? #
##########################################################

# This might be because there is a relationship between low MAF and high DAF
# Test whether there is a bias towards minor alleles with age having high variability 

# We need to now polarise with respect to minor allele
View(MAF_DAF)

# Create one file indicating the minor allele in each condition
Ctrl_minor_allele_table <- data.frame(SNP=rownames(MAF_DAF),Minor=MAF_DAF$Ctrl_minor_allele)
HS_minor_allele_table <- data.frame(SNP=rownames(MAF_DAF),Minor=MAF_DAF$HS_minor_allele)

# Save as tables
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Minor_Allele_Increased_Fraction/")
write.table(Ctrl_minor_allele_table,'Ctrl_minor_allele_table.txt',row.names = F,col.names = T,quote = F)
write.table(HS_minor_allele_table,'HS_minor_allele_table.txt',row.names = F,col.names = T,quote = F)

# The next script tests whether the median and skewness from QTL is different from non-QTL
# Tests whether the pattern is reproducible under subsets of the data

#################################################
##### End) Save working environment #############
#################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
save.image(file='6_Investigating_eQTL_and_veQTL_Allele_Age2.RData')



################# Old: non-downsampling method for veQTL #####################

# First let us reduce down the FDI lists to those with more than 1 gene 
# And sort these by MAF
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Derived_Allele_Increased_Fraction_Null")

# MAF Criteria.txt
Ctrl_FDI_SNP_MAF <- MAF[Ctrl_trans_per_veQTL_FDI$SNP,'Ctrl_MAF']
Ctrl_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  Ctrl_MAF_Criteria[quintile,'Min']=quantile(Ctrl_FDI_SNP_MAF,min)
  Ctrl_MAF_Criteria[quintile,'Max']=quantile(Ctrl_FDI_SNP_MAF,max)
  Ctrl_MAF_Criteria[quintile,'NumSNPstoSample']=round(length(Ctrl_FDI_SNP_MAF)/5)
}
write.table(Ctrl_MAF_Criteria,'Ctrl_trans-veQTL_MAF_Criteria.txt',row.names = F,col.names = T,quote = F)

HS_FDI_SNP_MAF <- MAF[HS_trans_per_veQTL_FDI$SNP,'HS_MAF']
HS_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  HS_MAF_Criteria[quintile,'Min']=quantile(HS_FDI_SNP_MAF,min)
  HS_MAF_Criteria[quintile,'Max']=quantile(HS_FDI_SNP_MAF,max)
  HS_MAF_Criteria[quintile,'NumSNPstoSample']=round(length(HS_FDI_SNP_MAF)/5)
}
write.table(HS_MAF_Criteria,'HS_trans-veQTL_MAF_Criteria.txt',row.names = F,col.names = T,quote = F)
View(HS_MAF_Criteria)

# non-veQTL_SNPs.txt - SNPs and their MAFs
# needs to have age and NOT be a veQTL
non_Ctrl_veQTL_aged<- intersect(setdiff(SNPs_position_map$SNP,Ctrl_trans_per_veQTL_FDI$SNP),unique(SNP_allele_age$SNP))
length(non_Ctrl_veQTL_aged)
# 121K of these 
non_Ctrl_trans_veQTL_by_MAF <- data.frame(SNP=non_Ctrl_veQTL_aged,MAF=MAF[non_Ctrl_veQTL_aged,'Ctrl_MAF'])
View(non_Ctrl_trans_veQTL_by_MAF)
write.table(non_Ctrl_trans_veQTL_by_MAF,'Ctrl_non-trans-veQTL_SNPs.txt',row.names = F,col.names = T,quote = F)
# Repeat for HS
non_HS_veQTL_aged<- intersect(setdiff(SNPs_position_map$SNP,HS_trans_per_veQTL_FDI$SNP),unique(SNP_allele_age$SNP))
length(non_HS_veQTL_aged)
# 95K of these 
non_HS_trans_veQTL_by_MAF <- data.frame(SNP=non_HS_veQTL_aged,MAF=MAF[non_HS_veQTL_aged,'HS_MAF'])
View(non_HS_trans_veQTL_by_MAF)
write.table(non_HS_trans_veQTL_by_MAF,'HS_non-trans-veQTL_SNPs.txt',row.names = F,col.names = T,quote = F)

# HS a little bit weird - check how many SNPs per MAF category can be selected for these SNPs
for (row in 1:nrow(HS_MAF_Criteria)){
  Min <-HS_MAF_Criteria[row,'Min']
  Max <-HS_MAF_Criteria[row,'Max']
  print(nrow(subset(non_HS_trans_veQTL_by_MAF, MAF>Min & MAF<Max)))
}
# 15K - 24K -> there should be a lot of variation to choose within the set
# It is not that the algorithm selects the same 5K SNPs every time...
# Though there is a more limited set of lowest MAFs than highest MAFs...
