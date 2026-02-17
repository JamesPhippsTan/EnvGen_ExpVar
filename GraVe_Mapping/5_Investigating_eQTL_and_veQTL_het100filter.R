# Investigating eQTL and veQTL
# Add a post-hoc filter for heterozygotes = 100

# Last Updated: 4/11/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls(all = T))

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
library(eulerr)
library(scales) 
library(ggh4x)
library(svglite)

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
#load(file='5_Investigating_eQTL_and_veQTL_het100.RData')

# Functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('gProfiler_Functions.R')
source('Quantile_Functions.R')
source('Variability_Functions.R')
source('QTL_Analysis_Functions.R')

#######################
##### Datasets ########
#######################

# Load veQTL result files
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
Ctrl_cis_veQTL <- read.table("Ctrl_cis_veQTL.txt", header = T,skipNul = T)
Ctrl_trans_veQTL <- read.table("Ctrl_trans_veQTL.0_0005.pval.txt", header = T,skipNul = T)
HS_cis_veQTL <- read.table("HS_cis_veQTL.txt", header = T,skipNul = T)
HS_trans_veQTL <- read.table("HS_trans_veQTL.0_0005.pval.txt", header = T,skipNul = T)

# Position file to re-map actual SNP names to trans-veQTL
SNPs_dummy_positions <- read.table("SNPs_dummy_positions.txt")
SNPs_original_positions <- read.table("SNPs_original_positions.txt")
SNPs_position_map <- merge(SNPs_dummy_positions[,2:3],SNPs_original_positions[,1:3],by='V3')
colnames(SNPs_position_map) <- c("SNP","POS","Actual_CHROM","Actual_POS")

# Gene body locations obtained from droseqtl.org - these match the dataset though are not necessarily the latest annotations on flybase
Gene_body_locations <-fread('droseQTLorg_genes_positions.csv',header = T,select = c(1,2,3,6),
                            col.names = c("GENE","GENE_START","GENE_END","GENE_CHR"))

# Load significant eQTL - all
Ctrl_cis_eqtl_sig <- read.table('Ctrl_cis_eqtl_sig.txt',header = T, fill=TRUE)
Ctrl_trans_eqtl_sig <- read.table('Ctrl_trans_eqtl_sig.txt',header = T, fill=TRUE)
HS_cis_eqtl_sig <- read.table('HS_cis_eqtl_sig.txt',header = T, fill=TRUE)
HS_trans_eqtl_sig <- read.table('HS_trans_eqtl_sig.txt',header = T, fill=TRUE)

# Significant eQTL by number of genes
Ctrl_cis_eqtl_number <- read.table('Ctrl_cis_eqtl_sig_number.txt',header = T, fill=TRUE)
Ctrl_trans_eqtl_number <- read.table('Ctrl_trans_eqtl_sig_number.txt',header = T, fill=TRUE)
HS_cis_eqtl_number <- read.table('HS_cis_eqtl_sig_number.txt',header = T, fill=TRUE)
HS_trans_eqtl_number <- read.table('HS_trans_eqtl_sig_number.txt',header = T, fill=TRUE)

# Load SNP stats
Ctrl_SNP_stats <- read.table("Dmel_Ctrl_final_geno_counts.frqx",header = T,sep='\t', check.names = FALSE)
Ctrl_SNP_stats <- add_MAF(Ctrl_SNP_stats)
HS_SNP_stats <- read.table("Dmel_HS_final_geno_counts.frqx",header = T,sep='\t', check.names = FALSE)
HS_SNP_stats <- add_MAF(HS_SNP_stats)
# Note: For most SNPs, A1 is the minor allele
MAF_initial <- data.frame(row.names = Ctrl_SNP_stats$SNP, 
                          Ctrl_MAF=Ctrl_SNP_stats$MAF,
                          HS_MAF=HS_SNP_stats$MAF,
                          Ctrl_minor_allele=Ctrl_SNP_stats$Minor_allele,
                          HS_minor_allele=HS_SNP_stats$Minor_allele,
                          MAF_Diff=(HS_SNP_stats$MAF-Ctrl_SNP_stats$MAF),
                          Avg_MAF =(Ctrl_SNP_stats$MAF+HS_SNP_stats$MAF)/2) 
MAF <- MAF_initial[SNPs_position_map$SNP,]

# Load SNPs with allele age - i.e., ancestral or derived
SNP_allele_age <- read.table("Snp_ancestral_state_results_onlyunambiguoushits.txt",header = T)
SNP_allele_age$SNP <- paste0(SNP_allele_age$arm,SNP_allele_age$snp_pos)
View(SNP_allele_age)
length(intersect(unique(SNP_allele_age$SNP),unique(SNPs_position_map$SNP)))
# 121683

# Whole-population mean and variability metrics to locate high-eQTL genes in the mean-variance plot
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
All_genes <- read.csv("Gene_MeanVar_Table.csv",row.names = 1)

#####################################
##### 1) Pre-processing #############
#####################################

# Remove cis-SNPs from trans-veQTL mapping results
Ctrl_trans_veQTL <- remove_cis_veQTL(trans_mapping_df = Ctrl_trans_veQTL, 
                                     SNP_positions = SNPs_position_map,
                                     gene_positions = Gene_body_locations,
                                     cis_window = 10000)
HS_trans_veQTL <- remove_cis_veQTL(trans_mapping_df = HS_trans_veQTL, 
                                   SNP_positions = SNPs_position_map,
                                   gene_positions = Gene_body_locations,
                                   cis_window = 10000)

# For cis, add SNP names (these names have been added to the trans results in the process of removing cis-veQTL)
Ctrl_cis_veQTL$SNP <- paste0(sub(23,"X",Ctrl_cis_veQTL$CHROM),Ctrl_cis_veQTL$POS)
HS_cis_veQTL$SNP <- paste0(sub(23,"X",HS_cis_veQTL$CHROM),HS_cis_veQTL$POS)

# For cis and trans, add gene-SNP pair names 
Ctrl_cis_veQTL$gene_SNP_pair <- paste0(Ctrl_cis_veQTL$GENE,"_",Ctrl_cis_veQTL$SNP)
HS_cis_veQTL$gene_SNP_pair <- paste0(HS_cis_veQTL$GENE,"_",HS_cis_veQTL$SNP)
Ctrl_trans_veQTL$gene_SNP_pair <- paste0(Ctrl_trans_veQTL$GENE,"_",Ctrl_trans_veQTL$SNP)
HS_trans_veQTL$gene_SNP_pair <- paste0(HS_trans_veQTL$GENE,"_",HS_trans_veQTL$SNP)

# Do the same for the eQTLs
Ctrl_cis_eqtl_sig$gene_SNP_pair <- paste0(Ctrl_cis_eqtl_sig$phenotype_id,"_",Ctrl_cis_eqtl_sig$variant_id)
HS_cis_eqtl_sig$gene_SNP_pair <- paste0(HS_cis_eqtl_sig$phenotype_id,"_",HS_cis_eqtl_sig$variant_id)
Ctrl_trans_eqtl_sig$gene_SNP_pair <- paste0(Ctrl_trans_eqtl_sig$phenotype_id,"_",Ctrl_trans_eqtl_sig$variant_id)
HS_trans_eqtl_sig$gene_SNP_pair <- paste0(HS_trans_eqtl_sig$phenotype_id,"_",HS_trans_eqtl_sig$variant_id)

######################################################################################
##### 2) Carry out multiple testing corrections to define significant veQTL ##########
######################################################################################

# Different p-value thresholds
p_val_cutoffs <- c(0.2,0.1,0.05)
for (p_val in p_val_cutoffs){
  print(p_val)
  
  # 1: BH correction for cis-QTL
  n_cis_tests <- nrow(Ctrl_cis_veQTL)
  n_cis_tests
  length(unique(Ctrl_cis_veQTL$GENE))
  length(unique(Ctrl_cis_veQTL$SNP))
  
  Ctrl_cis_veQTL$vfdr <- p.adjust(Ctrl_cis_veQTL$P,method = 'BH',n=n_cis_tests)
  Ctrl_cis_veQTL_sig <- subset(Ctrl_cis_veQTL,vfdr<p_val) 
  max(Ctrl_cis_veQTL_sig$P) 
  Ctrl_cis_veQTL_FDR0.05 <- max(subset(Ctrl_cis_veQTL,vfdr<0.05)$P) 
  Ctrl_cis_veQTL_FDR0.1 <- max(subset(Ctrl_cis_veQTL,vfdr<0.1)$P) 
  
  HS_cis_veQTL$vfdr <- p.adjust(HS_cis_veQTL$P,method = 'BH',n=n_cis_tests)
  HS_cis_veQTL_sig <- subset(HS_cis_veQTL,vfdr<p_val) 
  max(HS_cis_veQTL_sig$P) 
  HS_cis_veQTL_FDR0.05 <- max(subset(HS_cis_veQTL,vfdr<0.05)$P) 
  HS_cis_veQTL_FDR0.1 <- max(subset(HS_cis_veQTL,vfdr<0.1)$P) 
  
  # 2: BH correction for trans-QTL
  n_trans_tests <- 8763*383710-n_cis_tests  
  n_trans_tests
  length(unique(Ctrl_trans_veQTL$SNP))
  
  Ctrl_trans_veQTL$vfdr <- p.adjust(Ctrl_trans_veQTL$P,method = 'BH',n=n_trans_tests)
  Ctrl_trans_veQTL_sig <- subset(Ctrl_trans_veQTL,vfdr<p_val) 
  max(Ctrl_trans_veQTL_sig$P) # Critical raw pval = 1.07675e-08
  Ctrl_trans_veQTL_FDR0.05 <- max(subset(Ctrl_trans_veQTL,vfdr<0.05)$P) 
  Ctrl_trans_veQTL_FDR0.1 <- max(subset(Ctrl_trans_veQTL,vfdr<0.1)$P) 
  
  HS_trans_veQTL$vfdr <- p.adjust(HS_trans_veQTL$P,method = 'BH',n=n_trans_tests)
  HS_trans_veQTL_sig <- subset(HS_trans_veQTL,vfdr<p_val) 
  max(HS_trans_veQTL_sig$P) # Critical raw pval = 7.90163e-06
  HS_trans_veQTL_FDR0.05 <- max(subset(HS_trans_veQTL,vfdr<0.05)$P) 
  HS_trans_veQTL_FDR0.1 <- max(subset(HS_trans_veQTL,vfdr<0.1)$P) 
}

################################################################
#######  Investigate veQTL and their relation to eQTL ##########
################################################################

# Add SNP stats and info
Ctrl_cis_eqtl_sig <- add_SNP_stats(Ctrl_cis_eqtl_sig,
                                               SNP_stats_df = Ctrl_SNP_stats,
                                               merge_column = 'variant_id') 
HS_cis_eqtl_sig <- add_SNP_stats(HS_cis_eqtl_sig,
                                             SNP_stats_df = HS_SNP_stats,
                                             merge_column = 'variant_id')
Ctrl_trans_eqtl_sig <- add_SNP_stats(Ctrl_trans_eqtl_sig,
                                                 SNP_stats_df = Ctrl_SNP_stats,
                                                 merge_column = 'variant_id') 
HS_trans_eqtl_sig <- add_SNP_stats(HS_trans_eqtl_sig,
                                               SNP_stats_df = HS_SNP_stats,
                                               merge_column = 'variant_id') 

Ctrl_cis_veQTL_sig <- add_SNP_stats(Ctrl_cis_veQTL_sig,
                                                SNP_stats_df = Ctrl_SNP_stats,
                                                merge_column = 'SNP') 
HS_cis_veQTL_sig <- add_SNP_stats(HS_cis_veQTL_sig,
                                              SNP_stats_df = HS_SNP_stats,
                                              merge_column = 'SNP')
Ctrl_trans_veQTL_sig <- add_SNP_stats(Ctrl_trans_veQTL_sig,
                                                  SNP_stats_df = Ctrl_SNP_stats,
                                                  merge_column = 'SNP') 
HS_trans_veQTL_sig <- add_SNP_stats(HS_trans_veQTL_sig,
                                                SNP_stats_df = HS_SNP_stats,
                                                merge_column = 'SNP')


###############################
# (0) Impose het = 100 filter #
###############################

Ctrl_cis_eqtl_sig <- Ctrl_cis_eqtl_sig %>%  filter(`C(HET)` >= 100)
HS_cis_eqtl_sig <- HS_cis_eqtl_sig %>%  filter(`C(HET)` >= 100)
Ctrl_trans_eqtl_sig <- Ctrl_trans_eqtl_sig %>%  filter(`C(HET)` >= 100)
HS_trans_eqtl_sig <- HS_trans_eqtl_sig %>%  filter(`C(HET)` >= 100)

Ctrl_cis_veQTL_sig <- Ctrl_cis_veQTL_sig %>%  filter(`C(HET)` >= 100)
HS_cis_veQTL_sig <- HS_cis_veQTL_sig %>%  filter(`C(HET)` >= 100)
Ctrl_trans_veQTL_sig <- Ctrl_trans_veQTL_sig %>%  filter(`C(HET)` >= 100)
HS_trans_veQTL_sig <- HS_trans_veQTL_sig %>%  filter(`C(HET)` >= 100)

####################################################################################
# (1) How many significant gene-SNP veQTL pairs across the transcriptome-genome ####
####################################################################################

p_val_text <- paste0(' p < ',p_val)
print(paste0('Ctrl_cis_veQTL_sig',p_val_text))
print(nrow(Ctrl_cis_veQTL_sig)) # 11 - all over 100 in the HET category
print(paste0('HS_cis_veQTL_sig',p_val_text))
print(nrow(HS_cis_veQTL_sig)) # 148; 118 if more than 100 in the HET
print(paste0('Ctrl_trans_veQTL_sig',p_val_text))
print(nrow(Ctrl_trans_veQTL_sig)) # 724; 488 if more than 100 in the HET
print(paste0('HS_trans_veQTL_sig',p_val_text))
print(nrow(HS_trans_veQTL_sig)) # 531266: 380850 if more than 100 in the HET
 
print(paste0('Ctrl_cis_eqtl_sig',p_val_text))
print(nrow(Ctrl_cis_eqtl_sig)) # 86250
print(paste0('HS_cis_eqtl_sig',p_val_text))
print(nrow(HS_cis_eqtl_sig)) # 60189
print(paste0('Ctrl_trans_eqtl_sig',p_val_text))
print(nrow(Ctrl_trans_eqtl_sig)) # 702776
print(paste0('HS_trans_eqtl_sig',p_val_text))
print(nrow(HS_trans_eqtl_sig)) # 383026


#########################################################################################
# (1b) Dividing these in control-only shared, HS-only shared, and shared - pval #########
#########################################################################################

# Initialize the output dataframe
row_labels <- c("cis_FDR<0.2", "cis_FDR<0.1", "cis_FDR<0.05", "trans_FDR<0.2", "trans_FDR<0.1", "trans_FDR<0.05","cisandtrans_FDR<0.2", "cisandtrans_FDR<0.1", "cisandtrans_FDR<0.05")

veQTL_GxE_table <- data.frame(
  SharedCutoff = row_labels,
  GxE_Ctrl_genes = integer(length(row_labels)),
  GxE_HS_genes = integer(length(row_labels)),
  Shared_genes = integer(length(row_labels)),
  GxE_Ctrl_SNPs = integer(length(row_labels)),
  GxE_HS_SNPs = integer(length(row_labels)),
  Shared_SNPs = integer(length(row_labels)),
  stringsAsFactors = FALSE
)

### --------- CIS: FDR thresholds ---------
FDR_cutoffs <- c(0.2, 0.1, 0.05)
for (FDR in FDR_cutoffs) {
  row_name <- paste0("cis_FDR<", FDR)
  
  Ctrl_cis_cutoff <- subset(Ctrl_cis_veQTL, vfdr < FDR)
  HS_cis_cutoff <- subset(HS_cis_veQTL, vfdr < FDR)
  
  Shared_genes <- length(unique(c(
    intersect(Ctrl_cis_veQTL_sig$GENE, HS_cis_cutoff$GENE),
    intersect(HS_cis_veQTL_sig$GENE, Ctrl_cis_cutoff$GENE)
  )))
  
  GxE_Ctrl_genes <- length(setdiff(Ctrl_cis_veQTL_sig$GENE, HS_cis_cutoff$GENE))
  GxE_HS_genes   <- length(setdiff(HS_cis_veQTL_sig$GENE, Ctrl_cis_cutoff$GENE))
  
  Shared_SNPs <- length(unique(c(
    intersect(Ctrl_cis_veQTL_sig$SNP, HS_cis_cutoff$SNP),
    intersect(HS_cis_veQTL_sig$SNP, Ctrl_cis_cutoff$SNP)
  )))
  
  GxE_Ctrl_SNPs <- length(setdiff(Ctrl_cis_veQTL_sig$SNP, HS_cis_cutoff$SNP))
  GxE_HS_SNPs   <- length(setdiff(HS_cis_veQTL_sig$SNP, Ctrl_cis_cutoff$SNP))
  
  veQTL_GxE_table[veQTL_GxE_table$SharedCutoff == row_name, 2:7] <- c(
    GxE_Ctrl_genes, GxE_HS_genes, Shared_genes,
    GxE_Ctrl_SNPs, GxE_HS_SNPs, Shared_SNPs
  )
}

### --------- (2) TRANS: FDR thresholds ---------
for (FDR in FDR_cutoffs) {
  row_name <- paste0("trans_FDR<", FDR)
  
  Ctrl_trans_cutoff <- subset(Ctrl_trans_veQTL, vfdr < FDR)
  HS_trans_cutoff   <- subset(HS_trans_veQTL, vfdr < FDR)
  
  Ctrl_trans_cutoff_gene <- Ctrl_trans_cutoff$GENE
  HS_trans_cutoff_gene   <- HS_trans_cutoff$GENE
  Ctrl_trans_cutoff_snp  <- Ctrl_trans_cutoff$SNP
  HS_trans_cutoff_snp    <- HS_trans_cutoff$SNP
  
  Shared_genes <- length(unique(c(
    intersect(Ctrl_trans_veQTL_sig$GENE, HS_trans_cutoff_gene),
    intersect(HS_trans_veQTL_sig$GENE, Ctrl_trans_cutoff_gene)
  )))
  GxE_Ctrl_genes <- length(setdiff(Ctrl_trans_veQTL_sig$GENE, HS_trans_cutoff_gene))
  GxE_HS_genes   <- length(setdiff(HS_trans_veQTL_sig$GENE, Ctrl_trans_cutoff_gene))
  
  Shared_SNPs <- length(unique(c(
    intersect(Ctrl_trans_veQTL_sig$SNP, HS_trans_cutoff_snp),
    intersect(HS_trans_veQTL_sig$SNP, Ctrl_trans_cutoff_snp)
  )))
  GxE_Ctrl_SNPs <- length(setdiff(Ctrl_trans_veQTL_sig$SNP, HS_trans_cutoff_snp))
  GxE_HS_SNPs   <- length(setdiff(HS_trans_veQTL_sig$SNP, Ctrl_trans_cutoff_snp))
  
  veQTL_GxE_table[veQTL_GxE_table$SharedCutoff == row_name, 2:7] <- c(
    GxE_Ctrl_genes, GxE_HS_genes, Shared_genes,
    GxE_Ctrl_SNPs, GxE_HS_SNPs, Shared_SNPs
  )
}

### --------- (3) CISTRANS: FDR thresholds ---------
for (FDR in FDR_cutoffs) {
  row_name <- paste0("cisandtrans_FDR<", FDR)
  
  # New cutoff definitions
  Ctrl_cis_cutoff <- subset(Ctrl_cis_veQTL, vfdr < FDR)
  HS_cis_cutoff <- subset(HS_cis_veQTL, vfdr < FDR)
  Ctrl_trans_cutoff <- subset(Ctrl_trans_veQTL, vfdr < FDR)
  HS_trans_cutoff   <- subset(HS_trans_veQTL, vfdr < FDR)
  Ctrl_trans_cutoff_gene <- Ctrl_trans_cutoff$GENE
  HS_trans_cutoff_gene   <- HS_trans_cutoff$GENE
  Ctrl_trans_cutoff_snp  <- Ctrl_trans_cutoff$SNP
  HS_trans_cutoff_snp    <- HS_trans_cutoff$SNP
  
  # Get the intersects
  Shared_genes <- length(unique(c(
    intersect(c(Ctrl_cis_veQTL_sig$GENE,Ctrl_trans_veQTL_sig$GENE), 
              c(HS_cis_cutoff$GENE,HS_trans_cutoff_gene)),
    intersect(c(HS_cis_veQTL_sig$GENE,HS_trans_veQTL_sig$GENE), 
              c(Ctrl_cis_cutoff$GENE,Ctrl_trans_cutoff_gene))
  )))
  
  GxE_Ctrl_genes <- length(setdiff(c(Ctrl_cis_veQTL_sig$GENE,Ctrl_trans_veQTL_sig$GENE), 
                                   c(HS_cis_cutoff$GENE,HS_trans_cutoff_gene)))
  GxE_HS_genes   <- length(setdiff(c(HS_cis_veQTL_sig$GENE,HS_trans_veQTL_sig$GENE), 
                                   c(Ctrl_cis_cutoff$GENE,Ctrl_trans_cutoff_gene)))
  
  Shared_SNPs <- length(unique(c(
    intersect(c(Ctrl_cis_veQTL_sig$SNP,Ctrl_trans_veQTL_sig$SNP), 
              c(HS_cis_cutoff$SNP,HS_trans_cutoff_snp)),
    intersect(c(HS_cis_veQTL_sig$SNP,HS_trans_veQTL_sig$SNP), 
              c(Ctrl_cis_cutoff$SNP,Ctrl_trans_cutoff_snp))
  )))
  
  GxE_Ctrl_SNPs <- length(setdiff(c(Ctrl_cis_veQTL_sig$SNP,Ctrl_trans_veQTL_sig$SNP), 
                                  c(HS_cis_cutoff$SNP,HS_trans_cutoff_snp)))
  GxE_HS_SNPs   <- length(setdiff(c(HS_cis_veQTL_sig$SNP,HS_trans_veQTL_sig$SNP), 
                                  c(Ctrl_cis_cutoff$SNP,Ctrl_trans_cutoff_snp)))
  
  veQTL_GxE_table[veQTL_GxE_table$SharedCutoff == row_name, 2:7] <- c(
    GxE_Ctrl_genes, GxE_HS_genes, Shared_genes,
    GxE_Ctrl_SNPs, GxE_HS_SNPs, Shared_SNPs
  )
}

# View and save final summary table for this
View(veQTL_GxE_table)
#setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\GxE")
#write.csv(veQTL_GxE_table,'veQTL_GxE_table.csv')

# Turn it into a percentage
veQTL_GxE_table_percent <- veQTL_GxE_table
veQTL_GxE_table_percent[1:3,2:4] <- veQTL_GxE_table[1:3,2:4]*100/(sum(veQTL_GxE_table[1,2:4]))
veQTL_GxE_table_percent[4:6,2:4] <- veQTL_GxE_table[4:6,2:4]*100/(sum(veQTL_GxE_table[5,2:4]))
veQTL_GxE_table_percent[7:9,2:4] <- veQTL_GxE_table[7:9,2:4]*100/(sum(veQTL_GxE_table[9,2:4]))
veQTL_GxE_table_percent[1:3,5:7] <- veQTL_GxE_table[1:3,5:7]*100/(sum(veQTL_GxE_table[1,5:7]))
veQTL_GxE_table_percent[4:6,5:7] <- veQTL_GxE_table[4:6,5:7]*100/(sum(veQTL_GxE_table[5,5:7]))
veQTL_GxE_table_percent[7:9,5:7] <- veQTL_GxE_table[7:9,5:7]*100/(sum(veQTL_GxE_table[9,5:7]))
veQTL_GxE_table_percent[,2:7] <- round(veQTL_GxE_table_percent[,2:7],digits = 2)
View(veQTL_GxE_table_percent)
#write.csv(veQTL_GxE_table_percent,'veQTL_GxE_table_percent.csv')

# Make a plot of the percentages
veQTL_GxE_table_percent_separated <- veQTL_GxE_table_percent %>%
  separate(SharedCutoff, into = c("cis_or_trans", "shared_cutoff"), sep = "_")
# Shared is the opposite direction actually, turn the signs around
veQTL_GxE_table_percent_separated$shared_cutoff <- 
  factor(veQTL_GxE_table_percent_separated$shared_cutoff,
         levels = c("FDR<0.05", "FDR<0.1", "FDR<0.2")
  )
# Define the columns to plot
columns_to_plot <- colnames(veQTL_GxE_table_percent_separated)[3:8]

# Loop over columns
for (col_name in columns_to_plot) {
  p <- ggplot(veQTL_GxE_table_percent_separated,
              aes(x = shared_cutoff, y = .data[[col_name]],
                  color = cis_or_trans, group = cis_or_trans)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    theme_classic() +
    labs(
      title = 'veQTL',
      x = "Significance cutoff in other condition",
      y = gsub('_',' ',paste0(col_name, " (%)")),
      color = "Regulation Type"
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  # Save plot
  ggsave(filename = paste0(col_name, "_veQTL_dotplot_het100.svg"),
         plot = p, width = 6, height = 4, dpi = 300)
}

# At Het = 100, still similar -> mostly HS-specific genes and SNPs


#########################
# (2) eGenes and vGenes #
#########################

##############################################################################
# (2a) How many vGenes and eGenes in the transcriptome? Are there universal vGenes or eGenes?
##############################################################################

cis_vGenes_Ctrl <- unique(Ctrl_cis_veQTL_sig$GENE)
cis_vGenes_HS <- unique(HS_cis_veQTL_sig$GENE)
trans_vGenes_Ctrl <- unique(Ctrl_trans_veQTL_sig$GENE)
trans_vGenes_HS <- unique(HS_trans_veQTL_sig$GENE)
cis_eGenes_Ctrl <- unique(Ctrl_cis_eqtl_sig$phenotype_id)
cis_eGenes_HS <- unique(HS_cis_eqtl_sig$phenotype_id)
trans_eGenes_Ctrl <- unique(Ctrl_trans_eqtl_sig$phenotype_id)
trans_eGenes_HS <- unique(HS_trans_eqtl_sig$phenotype_id)

eGenes_Ctrl <- c(cis_eGenes_Ctrl,trans_eGenes_Ctrl)
eGenes_HS <- c(cis_eGenes_HS,trans_eGenes_HS)
vGenes_Ctrl <- c(cis_vGenes_Ctrl,trans_vGenes_Ctrl)
vGenes_HS <- c(cis_vGenes_HS,trans_vGenes_HS)

# Overlap plots, i.e., horizontal barplots 
# Create separate data frames for each category with pre-calculated overlaps
# FDR>0.05 in one condition, FDR>0.1 in the other
# Combined cis and trans
# Done for all eQTL and veQTL
data_eQTL <- data.frame(
  Partition = c("No", "Ctrl only", "Ctrl and HS", "HS only"),
  Value = c(8763-571-7531-50, 571, 7531, 50) 
) %>%
  mutate(
    Percent = Value / sum(Value),
    Partition = factor(Partition, levels = c("Ctrl and HS", "HS only", "Ctrl only", "No"))
  )

data_veQTL <- data.frame(
  Partition = c("No", "Ctrl only", "Ctrl and HS", "HS only"),
  Value = c(8763-0-1230-6824, 0, 1230, 6824) 
) %>%
  mutate(
    Percent = Value / sum(Value),
    Partition = factor(Partition, levels = c("Ctrl and HS", "HS only", "Ctrl only", "No"))
  )

# Colours
veQTL_colors <- c(
  "HS only" = "#9651EC",
  "No veQTL" = "#2D0B58",
  "Ctrl and HS" = "#CFB2F3",
  "Ctrl only" = "#611BB8"
)
eQTL_colors <- c(
  "HS only" = "#4F8E4D",
  "No eQTL" = "#192F19",
  "Ctrl and HS" = "#77BB75",
  "Ctrl only" = "#2B5D29"
)

# Create overlap plots
eQTL_plot <- create_overlap_plot(data_eQTL,eQTL_colors,'Number of genes with eQTL')
eQTL_plot
veQTL_plot <- create_overlap_plot(data_veQTL,veQTL_colors,'Number of genes with veQTL')
veQTL_plot

# Save plots
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Overlap_Ctrl_HS_eQTL_veQTL")
# Number of eQTL and veQTL per gene (with at least 1 eQTL or veQTL)
#ggsave(plot = eQTL_plot, filename='eGene_overlap_plot.svg', width = 3, height = 1.5, dpi = 300)
ggsave(plot = veQTL_plot, filename='vGene_overlap_plot_het100.svg', width = 3, height =1.5, dpi = 300)

#########################################################
# (2b) How many SNPs per Gene? Are there hotspot genes? #
#########################################################

# eGenes 
eQTL_hotspot_genes <- merge(Ctrl_cis_eqtl_number,HS_cis_eqtl_number,all = T,by='phenotype_id')
eQTL_hotspot_genes <- merge(eQTL_hotspot_genes,Ctrl_trans_eqtl_number,all = T,by='phenotype_id')
eQTL_hotspot_genes <- merge(eQTL_hotspot_genes,HS_trans_eqtl_number,all = T,by='phenotype_id')
colnames(eQTL_hotspot_genes) <- c('gene_id','Ctrl_cis_eQTL_#','HS_cis_eQTL_#','Ctrl_trans_eQTL_#','HS_trans_eQTL_#')

# vGenes 
veQTL_hotspot_genes <- merge(veQTL_number_per_gene(Ctrl_cis_veQTL_sig),veQTL_number_per_gene(HS_cis_veQTL_sig),all = T,by='GENE')
veQTL_hotspot_genes <- merge(veQTL_hotspot_genes,veQTL_number_per_gene(Ctrl_trans_veQTL_sig),all = T,by='GENE')
veQTL_hotspot_genes <- merge(veQTL_hotspot_genes,veQTL_number_per_gene(HS_trans_veQTL_sig),all = T,by='GENE')
colnames(veQTL_hotspot_genes) <- c('gene_id','Ctrl_cis_veQTL_#','HS_cis_veQTL_#','Ctrl_trans_veQTL_#','HS_trans_veQTL_#')

# Get number of each type of QTL per gene
QTL_number_per_gene <- merge(eQTL_hotspot_genes,veQTL_hotspot_genes,by = 'gene_id',all = T)
View(QTL_number_per_gene)
# Include genes without any QTL as well
Genes_with_no_QTL <- setdiff(rownames(All_genes),QTL_number_per_gene$gene_id)
# Create a new tibble with 0s for other columns
new_rows <- tibble(gene_id = Genes_with_no_QTL) %>%
  mutate(across(-gene_id, ~ 0))
# Combine the original df with the new rows
QTL_number_per_gene_final <- bind_rows(QTL_number_per_gene, new_rows)
# Turn NAs into 0s
QTL_number_per_gene_final[is.na(QTL_number_per_gene_final)] <- 0 


# Summary tables - includes the total number of genes with at least 1 eQTL or veQTL
eQTL_hotspot_genes_summary <- as.data.frame(t(sapply(eQTL_hotspot_genes[,2:ncol(eQTL_hotspot_genes)], summary_with_count)))
eQTL_hotspot_genes_summary
veQTL_hotspot_genes_summary <- as.data.frame(t(sapply(veQTL_hotspot_genes[,2:ncol(veQTL_hotspot_genes)], summary_with_count)))
veQTL_hotspot_genes_summary
QTL_number_per_gene_summary <- rbind(eQTL_hotspot_genes_summary,veQTL_hotspot_genes_summary)
View(QTL_number_per_gene_summary)

# Comparing with Huang et al., 2015, there is a range of 1-7 independent eQTL and veQTL per gene
# The median numbers fall within this range though the maximum numbers far surpass this


###############################################################################################
##### (2e) eQTL and veQTL number vs variability and mean rank in the whole transcriptome ######
###############################################################################################

All_genes <- rownames_to_column(All_genes,'gene_id')
hotspot_genes_All_genes <- merge(QTL_number_per_gene_final,All_genes,all=T,by='gene_id')
View(hotspot_genes_All_genes)

# Make plots
for (QTL_type in c("cis_eQTL_#","cis_veQTL_#","trans_eQTL_#","trans_veQTL_#")){
  mean_plot_list <- list()
  var_plot_list <- list()
  for (condition in c("Ctrl","HS")){
    
    # Last, loop for each condition
    var_category <- paste0(condition,"_VST_MAD")
    mean_category <- paste0(condition,"_VST_Mean")
    QTL_category <- paste0(condition,"_",QTL_type)
    
    # Omit 0s
    df_subset <- hotspot_genes_All_genes %>%
      dplyr::select(!!sym(var_category), !!sym(mean_category), !!sym(QTL_category)) %>%
      dplyr::filter(if_all(everything(), ~ !is.na(.) & . != 0))
    
    # Do not omit 0s
    df_subset <- hotspot_genes_All_genes %>%
      dplyr::select(!!sym(var_category), !!sym(mean_category), !!sym(QTL_category))
    
    df_subset[[var_category]] <- as.numeric(df_subset[[var_category]])
    df_subset[[mean_category]] <- as.numeric(df_subset[[mean_category]])
    df_subset[[QTL_category]] <- as.numeric(df_subset[[QTL_category]])
    
    # Extract Spearman correlation coefficient and p-value
    spearman_test <- cor.test(y=df_subset[[var_category]], x=df_subset[[QTL_category]], method = "spearman")
    spearman_coef <- spearman_test$estimate
    spearman_p_value <- spearman_test$p.value
    spearman_test_result <- paste0("Spearman rho = ",format_statistic(spearman_coef),
                                   " p = ",format_statistic(spearman_p_value),"\n",
                                   "Number of genes = ",nrow(df_subset))
    
    # Repeat spearman for mean
    mean_spearman_test <- cor.test(y=df_subset[[mean_category]], x=df_subset[[QTL_category]], method = "spearman")
    mean_spearman_coef <- mean_spearman_test$estimate
    mean_spearman_p_value <- mean_spearman_test$p.value
    mean_spearman_test_result <-paste0("Spearman rho = ",format_statistic(mean_spearman_coef),
                                       " p = ",format_statistic(mean_spearman_p_value),"\n",
                                       "Number of genes = ",nrow(df_subset))
    
    # Prepare label text
    Var_label <- paste0(condition," transcript level variability")
    Mean_label <- paste0(condition," mean transcript level")
    QTL_cat_cleaned <- gsub('_',' ',QTL_category)
    QTL_cat_cleaned <- gsub('trans ','trans-',QTL_cat_cleaned)
    QTL_cat_cleaned <- gsub('cis ','cis-',QTL_cat_cleaned)
    x_label <- paste0('log (',QTL_cat_cleaned,' + 1)')
    
    # Plot
    var_plot <- ggplot(df_subset, aes(y = !!sym(var_category), x = log(!!sym(QTL_category)+1))) +
      geom_point(alpha=0.5,color='darkgrey') +
      geom_smooth(method = "lm", col = "blue",linetype='dashed')+
      annotate("text", y = max(df_subset[[var_category]]), x = max(log(1+df_subset[[QTL_category]])), label = spearman_test_result, size = 3, hjust = 1,vjust = 1)+
      ylab(Var_label)+
      xlab(x_label)+
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            plot.title = element_text(size = 7.5))+
      theme_classic()
    
    mean_plot <- ggplot(df_subset, aes(y = !!sym(mean_category), x = log(!!sym(QTL_category)+1))) +
      geom_point(alpha=0.5,color='darkgrey') +
      geom_smooth(method = "lm", col = "blue",linetype='dashed')+
      annotate("text", y = max(df_subset[[mean_category]]), x = max(log(1+df_subset[[QTL_category]])), label = mean_spearman_test_result, size = 3, hjust = 1,vjust = 1)+
      ylab(Mean_label)+
      xlab(x_label)+
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            plot.title = element_text(size = 7.5))+
      theme_classic()
    
    mean_plot_list[[mean_category]] <- mean_plot
    var_plot_list[[var_category]] <- var_plot
  }
  # Save the plots
  setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\eQTL_veQTL_meanvar_plots")
  mean_plots <- wrap_plots(mean_plot_list, nrow = 2)  # Arrange in 2 columns
  print(mean_plots)
  #ggsave(plot = mean_plots,filename = paste0("MeanCor_",QTL_category,"plot.svg"),width=3,height=6,dpi=300)
  var_plots <- wrap_plots(var_plot_list, nrow = 2)  # Arrange in 2 columns
  print(var_plots)
  #ggsave(plot = var_plots,filename = paste0("VarCor_",QTL_category,"plot.svg"),width=3,height=6,dpi=300)
}

# cis-eQTL number
# Voom - Small (in HS, insignificant) positive with variability and mean
# VST - Larger positive with variability and small positive with mean

# cis-veQTL number
# Voom - Ctrl-cis: 3 high-variability genes
# VST - Ctrl-cis: 2 high-variability and one super low-variability

# trans-eQTL number
# Voom - Small negative with variability and positive with mean
# VST - Larger positive with variability and negative with mean

# trans-veQTL number (comment on Ctrl_trans)
# Voom - Positive with variability and negative with mean
# VST - Weaker positive with variability and negative with mean

# The cis- and trans- eQTL plots provide evidence that at least some of the high-mean high-variance genes in the VST plot are those with higher-than-normal numbers of eQTL
# The trans-veQTL plots provide evidence that the number of trans-veQTL scales with the raw variance

#######################
# (3) eQTLs and vQTLs #
#######################

#######################################################
# (3a) Are there SNPs shared across the conditions? ###
#######################################################

# eQTLs
cis_eQTL_Ctrl <- unique(Ctrl_cis_eqtl_sig$SNP)
trans_eQTL_Ctrl <- unique(Ctrl_trans_eqtl_sig$SNP)
cis_eQTL_HS <- unique(HS_cis_eqtl_sig$SNP)
trans_eQTL_HS <- unique(HS_trans_eqtl_sig$SNP)

# vQTLs
cis_veQTL_Ctrl <- unique(Ctrl_cis_veQTL_sig$SNP)
trans_veQTL_Ctrl <- unique(Ctrl_trans_veQTL_sig$SNP)
cis_veQTL_HS <- unique(HS_cis_veQTL_sig$SNP)
trans_veQTL_HS <- unique(HS_trans_veQTL_sig$SNP)

eQTL_Ctrl <- c(cis_eQTL_Ctrl,trans_eQTL_Ctrl)
eQTL_HS <- c(cis_eQTL_HS,trans_eQTL_HS)
veQTL_Ctrl <- c(cis_veQTL_Ctrl,trans_veQTL_Ctrl)
veQTL_HS <- c(cis_veQTL_HS,trans_veQTL_HS)

# Horizontal barplot visualisation
data_eQTL_SNPs <- data.frame(
  Partition = c("Non-eQTL", "Ctrl-only eQTL", "Ctrl and HS eQTL", "HS-only eQTL"),
  Value = c(383.710-62.929-8.968-237.976, 62.929, 237.976,8.968) 
) %>%
  mutate(
    Percent = Value / sum(Value),
    Partition = factor(Partition, levels = c("Ctrl and HS eQTL", "HS-only eQTL", "Ctrl-only eQTL", "Non-eQTL"))
  )

data_veQTL_SNPs <- data.frame(
  Partition = c("Non-veQTL", "Ctrl-only veQTL", "Ctrl and HS veQTL", "HS-only veQTL"),
  Value = c(383.710-0.063-1.848-88.119, 0.063, 1.848,88.119) 
) %>%
  mutate(
    Percent = Value / sum(Value),
    Partition = factor(Partition, levels = c("Ctrl and HS veQTL", "HS-only veQTL", "Ctrl-only veQTL", "Non-veQTL"))
  )

# Colours
eQTL_colors <- c(
  "HS-only eQTL" = "#4F8E4D",
  "Non-eQTL" = "#192F19",
  "Ctrl and HS eQTL" = "#77BB75",
  "Ctrl-only eQTL" = "#2B5D29"
)
veQTL_colors <- c(
  "HS-only veQTL" = "#9651EC",
  "Non-veQTL" = "#2D0B58",
  "Ctrl and HS veQTL" = "#CFB2F3",
  "Ctrl-only veQTL" = "#611BB8"
)

# Create overlap plots
eQTL_SNP_plot <- create_overlap_plot(data_eQTL_SNPs,eQTL_colors,
                                     expression(paste("Number of SNPs / x", 10^3)))
eQTL_SNP_plot
veQTL_SNP_plot <- create_overlap_plot(data_veQTL_SNPs,veQTL_colors,
                                      expression(paste("Number of SNPs / x", 10^3)))
veQTL_SNP_plot

# Save plots
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Overlap_Ctrl_HS_eQTL_veQTL")
# Number of eQTL and veQTL per gene (with at least 1 eQTL or veQTL)
#ggsave(plot = eQTL_SNP_plot, filename='eQTL_SNP_overlap_plot.svg', width = 3, height = 1.5, dpi = 300)
ggsave(plot = veQTL_SNP_plot, filename='veQTL_SNP_overlap_plot_het100.svg', width = 3, height =1.5, dpi = 300)

# Eliminates 5.5% of veQTL

##############################################################################
# (3b) How many genes per SNP? Are there hotspot SNPs that regulate many genes?
##############################################################################

# Get the number of hotspot QTL for each test within one table
hotspot_SNPs <- merge(data.frame(table(Ctrl_cis_eqtl_sig$SNP)),data.frame(table(Ctrl_trans_eqtl_sig$SNP)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(HS_cis_eqtl_sig$SNP)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(HS_trans_eqtl_sig$SNP)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(Ctrl_cis_veQTL_sig$SNP)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(HS_cis_veQTL_sig$SNP)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(Ctrl_trans_veQTL_sig$SNP)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(HS_trans_veQTL_sig$SNP)),all = T,by='Var1')
colnames(hotspot_SNPs) <- c('SNP','Ctrl_cis_eQTLs','Ctrl_trans_eQTLs','HS_cis_eQTLs','HS_trans_eQTLs','Ctrl_cis_veQTLs','HS_cis_veQTLs','Ctrl_trans_veQTLs','HS_trans_veQTLs')
hotspot_SNPs_summary <- as.data.frame(t(sapply(hotspot_SNPs[,2:ncol(hotspot_SNPs)], summary_with_count)))
View(hotspot_SNPs_summary)

# What determines the number of genes per SNP? 
HS_trans_veQTL_SNP_stats <- merge(hotspot_SNPs[,c('SNP','HS_trans_veQTLs')],HS_SNP_stats,by = 'SNP',all.x = )
# MAF?
ggplot(data=HS_trans_veQTL_SNP_stats,aes(x=MAF,y=log(HS_trans_veQTLs))) + geom_point()
cor.test(HS_trans_veQTL_SNP_stats$MAF,HS_trans_veQTL_SNP_stats$HS_trans_veQTLs,method='spearman')
# Tiny negative correlation
ggplot(data=HS_trans_veQTL_SNP_stats,aes(x=`C(HET)`,y=log(HS_trans_veQTLs))) + geom_point()
cor.test(HS_trans_veQTL_SNP_stats$`C(HET)`,HS_trans_veQTL_SNP_stats$HS_trans_veQTLs,method='spearman')
# Larger negative correlation
# How about the top 55...
summary(subset(HS_trans_veQTL_SNP_stats,HS_trans_veQTLs>400)$`C(HET)`)
# Median 125 at 200
# Median 101 at 300
# Median max 90 at 400, meaning the top 4 veQTL are low-confidence...
# The 1011 one is a low confidence-SNP
summary(subset(HS_trans_veQTL_SNP_stats,HS_trans_veQTLs>200)$`C(HET)`)


# There have been no comparison of the number of genes sharing eQTL/veQTL in Huang et al., 2015
# Save for the number of genetically-correlated transcripts based on MMC but not on eQTL/veQTL number
# Thus, these are new results 
# The strongest HS hotspot veQTL affects more genes than the strongest Ctrl veQTL
# Plot the pleiotropy per eQTL and veQTL category 
for (SNP_column in c('Ctrl_cis_eQTLs','Ctrl_trans_eQTLs','HS_cis_eQTLs','HS_trans_eQTLs')){
  data <- data.frame(value = hotspot_SNPs[, SNP_column])
  data <- na.omit(data)
  binwidth <- ifelse(max(data)<10,1,
                     (max(data$value) - min(data$value)) / max(data$value))  
  eQTL_pleiotropy <- ggplot(data) +
    geom_histogram(
      aes(x = value),
      binwidth = binwidth,
      center = 0,
      bins = 10,
      position = "identity",
      color = "black",
      fill = "#4F8E4D"
    ) +
    theme_classic() +
    labs(
      title = gsub('s e','s-e',gsub('_',' ', SNP_column)),
      subtitle = paste0(nrow(data),' SNPs'),
      x = "# of regulated genes per SNP",
      y = '# of SNPs (log scale)'
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )+
    scale_y_continuous(
      trans = pseudo_log_trans(base = 10),
      breaks = c(0, 1, 10, 100, 1000, 10000)
    )
  setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs")
  #ggsave(plot = eQTL_pleiotropy,filename = paste0("SNP_pleiotropy,",SNP_column,"_plot.svg"),height=3,width=3,dpi=300)
}

# Plot the pleiotropy per eQTL and veQTL category 
for (SNP_column in c('Ctrl_cis_veQTLs','Ctrl_trans_veQTLs','HS_cis_veQTLs','HS_trans_veQTLs')){
  data <- data.frame(value = hotspot_SNPs[, SNP_column])
  data <- na.omit(data)
  binwidth <- ifelse(max(data)<10,1,
                     (max(data$value) - min(data$value)) / max(data$value))
  veQTL_pleiotropy <- ggplot(data) +
    geom_histogram(
      aes(x = value),
      bins = 10,
      binwidth = binwidth,
      center = 0,
      position = "identity",
      color = "black",
      fill= '#611BB8'
    ) +
    theme_classic() +
    labs(
      title = gsub('s ve','s-ve',gsub('_',' ', SNP_column)),
      subtitle = paste0(nrow(data),' SNPs'),
      x = "# of regulated genes per SNP",
      y = "# of SNPs (log scale)"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )+
    scale_y_continuous(
      trans = pseudo_log_trans(base = 10),
      breaks = c(0, 1, 10, 100, 1000, 10000)
    )
  setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs")
  #ggsave(plot = veQTL_pleiotropy,filename = paste0("SNP_pleiotropy,",SNP_column,"_plot.svg"),height=3,width=3,dpi=300)
}

###############################################################################
# (3c) Do the SNP sets have biases in minor allele frequency? eQTL vs veQTL ###
###############################################################################

# eQTL list
eQTL_list <- list(
  All = unique(c(Ctrl_cis_eqtl_sig$SNP,
                 Ctrl_trans_eqtl_sig$SNP,
                 HS_cis_eqtl_sig$SNP,
                 HS_trans_eqtl_sig$SNP)),
  All_cis = unique(c(Ctrl_cis_eqtl_sig$SNP,
                     HS_cis_eqtl_sig$SNP)),
  All_trans = unique(c(Ctrl_trans_eqtl_sig$SNP,
                       HS_trans_eqtl_sig$SNP))
)

# veQTL list
veQTL_list <- list(
  All = unique(c(Ctrl_cis_veQTL_sig$SNP,
                 Ctrl_trans_veQTL_sig$SNP,
                 HS_cis_veQTL_sig$SNP,
                 HS_trans_veQTL_sig$SNP)),
  All_cis = unique(c(Ctrl_cis_veQTL_sig$SNP,
                     HS_cis_veQTL_sig$SNP)),
  All_trans = unique(c(Ctrl_trans_veQTL_sig$SNP,
                       HS_trans_veQTL_sig$SNP))
)

# Plot eQTL and veQTL MAF histograms next to each other
compare_MAFs_eQTL_veQTL_het <- function(eQTL_datasets, veQTL_datasets, MAF_df) {
  # MAF column selector
  get_maf_column <- function(name) {
    if (grepl("^Ctrl", name)) return("Ctrl_MAF")
    if (grepl("^HS", name)) return("HS_MAF")
    return("Avg_MAF")  # For "All" sets
  }
  
  results <- data.frame(
    Comparison = character(),
    KS_p = numeric(),
    T_p = numeric(),
    medianDiff = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(eQTL_datasets)) {
    if (!name %in% names(veQTL_datasets)) next
    
    message("Processing: ", name)
    maf_col <- get_maf_column(name)
    
    eqtl_snps <- unique(eQTL_datasets[[name]])
    veqtl_raw <- veQTL_datasets[[name]]
    veqtl_snps <- if (is.data.frame(veqtl_raw)) unique(veqtl_raw$SNP) else unique(veqtl_raw)
    
    # Filter to SNPs present in MAF_df
    eqtl_snps <- eqtl_snps[eqtl_snps %in% rownames(MAF_df)]
    veqtl_snps <- veqtl_snps[veqtl_snps %in% rownames(MAF_df)]
    
    eqtl_maf <- MAF_df[eqtl_snps, maf_col, drop = TRUE]
    veqtl_maf <- MAF_df[veqtl_snps, maf_col, drop = TRUE]
    
    # Remove NA values
    eqtl_maf <- eqtl_maf[!is.na(eqtl_maf)]
    veqtl_maf <- veqtl_maf[!is.na(veqtl_maf)]
    
    if (length(eqtl_maf) < 3 || length(veqtl_maf) < 3) next
    
    # Statistical tests
    ks_res <- ks.test(eqtl_maf, veqtl_maf)
    w_res <- wilcox.test(eqtl_maf, veqtl_maf)
    median_eqtl <- median(eqtl_maf)
    median_veqtl <- median(veqtl_maf)
    
    # Format p-values
    format_P_handle0 <- function(p) {
      if (is.na(p)) return("NA")
      else if (p == 0) return("2.2e-16")
      else return(formatC(p, format = "e", digits = 2))
    }
    ks_p <- format_P_handle0(ks_res$p.value)
    wilcox_p <- format_P_handle0(w_res$p.value)
    
    results <- rbind(results, data.frame(
      Comparison = name,
      KS_p = ks_p,
      Wilcox_p = wilcox_p,
      eQTL_median_MAF = median_eqtl,
      veQTL_median_MAF = median_veqtl,
      eQTL_n = length(eqtl_maf),
      veQTL_n = length(veqtl_maf)
    ))
    
    # Prepare plot data
    plot_data <- rbind(
      data.frame(Group = "eQTL", MAF = eqtl_maf),
      data.frame(Group = "veQTL", MAF = veqtl_maf)
    )
    
    # Bin and normalize within each group
    bin_breaks <- seq(0, 0.5, length.out = 11)
    plot_data_binned <- plot_data %>%
      mutate(Bin = cut(MAF, breaks = bin_breaks, include.lowest = TRUE, right = FALSE)) %>%
      group_by(Group, Bin) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Group) %>%
      mutate(Percent = Count / sum(Count) * 100)
    
    # Plot with all requested changes
    p_legend <- ggplot(plot_data_binned, aes(x = Bin, y = Percent, fill = Group)) +
      geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
      theme_classic() +
      labs(
        title = paste0("Kolmogorov-Smirnoff p = ", ks_p, " | Wilcox p = ", wilcox_p),
        x = "MAF",
        y = "Percent of all QTL"
      ) +
      scale_fill_manual(
        values = c("eQTL" = "#4F8E4D", "veQTL" = "#611BB8"),  # Updated eQTL color
        labels = c(
          paste0("eQTL (n = ", length(eqtl_maf), ", median = ", round(median_eqtl, 3), ")"),
          paste0("veQTL (n = ", length(veqtl_maf), ", median = ", round(median_veqtl, 3), ")")
        ),
        name = NULL) + theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = c(1, 1),  # Legend inside plot (top-right)
        legend.justification = c(1, 1),  # Anchor to top-right corner
        legend.background = element_rect(fill = "white", color = "black"),  # Legible background
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
      )
    p <- p_legend + theme(legend.position = "none")
    setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\MAF")
    ggsave(filename = paste0(name, "_MAF_comparison_het100.svg"), plot = p, width = 4.7, height = 4,dpi=300)
    ggsave(filename = paste0(name, "_MAF_comparison_het100legend.svg"), plot = p_legend, width = 4.7, height = 4,dpi=300)
  }
  return(results)
}

MAFs_eQTL_veQTL <- compare_MAFs_eQTL_veQTL_het(eQTL_datasets=eQTL_list, 
                                           veQTL_datasets=veQTL_list, 
                                           MAF_df=MAF)

# At het = 100, the tail is erased but the median is still lower
# Still different distributions with a lean towards low-frequency