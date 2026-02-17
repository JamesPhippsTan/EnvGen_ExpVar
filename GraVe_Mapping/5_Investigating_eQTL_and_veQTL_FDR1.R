# Investigating eQTL and veQTL
# 0.01 FDR as the cutoff

# Last Updated: 22/10/25

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
#load(file='5_Investigating_eQTL_and_veQTL_FDR1.RData')

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
p_val_cutoffs <- c(0.1,0.05,0.01)
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
  Ctrl_cis_veQTL_FDR0.01 <- max(subset(Ctrl_cis_veQTL,vfdr<0.01)$P) 
  
  HS_cis_veQTL$vfdr <- p.adjust(HS_cis_veQTL$P,method = 'BH',n=n_cis_tests)
  HS_cis_veQTL_sig <- subset(HS_cis_veQTL,vfdr<p_val) 
  max(HS_cis_veQTL_sig$P) 
  HS_cis_veQTL_FDR0.05 <- max(subset(HS_cis_veQTL,vfdr<0.05)$P) 
  HS_cis_veQTL_FDR0.01 <- max(subset(HS_cis_veQTL,vfdr<0.01)$P) 
  
  # 2: BH correction for trans-QTL
  n_trans_tests <- 8763*383710-n_cis_tests  
  n_trans_tests
  length(unique(Ctrl_trans_veQTL$SNP))
  
  Ctrl_trans_veQTL$vfdr <- p.adjust(Ctrl_trans_veQTL$P,method = 'BH',n=n_trans_tests)
  Ctrl_trans_veQTL_sig <- subset(Ctrl_trans_veQTL,vfdr<p_val) 
  max(Ctrl_trans_veQTL_sig$P) # Critical raw pval = 1.07675e-08
  Ctrl_trans_veQTL_FDR0.05 <- max(subset(Ctrl_trans_veQTL,vfdr<0.05)$P) 
  Ctrl_trans_veQTL_FDR0.01 <- max(subset(Ctrl_trans_veQTL,vfdr<0.01)$P) 
  
  HS_trans_veQTL$vfdr <- p.adjust(HS_trans_veQTL$P,method = 'BH',n=n_trans_tests)
  HS_trans_veQTL_sig <- subset(HS_trans_veQTL,vfdr<p_val) 
  max(HS_trans_veQTL_sig$P) # Critical raw pval = 7.90163e-06
  HS_trans_veQTL_FDR0.05 <- max(subset(HS_trans_veQTL,vfdr<0.05)$P) 
  HS_trans_veQTL_FDR0.01 <- max(subset(HS_trans_veQTL,vfdr<0.01)$P) 
  
  # Effect size (correlation slope) distributions
  cor_plot_list <- list()
  
  cor_plot_list[[1]] <-ggplot(Ctrl_cis_veQTL_sig,aes(y=abs(COR),x=P))+geom_point()+theme_classic()+geom_vline(xintercept = Ctrl_cis_veQTL_FDR0.05, color = "red", linetype = "dashed")+
    geom_vline(xintercept = Ctrl_cis_veQTL_FDR0.01, color = "red", linetype = "dashed")+
    ggtitle('Ctrl_cis_veQTL_sig_FDR1')
  cor_plot_list[[2]] <-ggplot(HS_cis_veQTL_sig,aes(y=abs(COR),x=P))+geom_point() +theme_classic()+geom_vline(xintercept = HS_cis_veQTL_FDR0.05, color = "red", linetype = "dashed")+
    geom_vline(xintercept = HS_cis_veQTL_FDR0.01, color = "red", linetype = "dashed")+
    ggtitle('HS_cis_veQTL_sig_FDR1')
  cor_plot_list[[3]] <-ggplot(Ctrl_trans_veQTL_sig,aes(y=abs(COR),x=P))+geom_point() +theme_classic()+geom_vline(xintercept = Ctrl_trans_veQTL_FDR0.05, color = "red", linetype = "dashed")+
    geom_vline(xintercept = Ctrl_trans_veQTL_FDR0.01, color = "red", linetype = "dashed")+
    ggtitle('Ctrl_trans_veQTL_sig_FDR1')
  cor_plot_list[[4]] <-ggplot(HS_trans_veQTL_sig,aes(y=abs(COR),x=P))+geom_point() +theme_classic()+geom_vline(xintercept = HS_trans_veQTL_FDR0.05, color = "red", linetype = "dashed")+
    geom_vline(xintercept = HS_trans_veQTL_FDR0.01, color = "red", linetype = "dashed")+
    ggtitle('HS_trans_veQTL_sig_FDR1')
  
  cor_plots <- cor_plot_list[[1]] + cor_plot_list[[2]] + cor_plot_list[[3]] + cor_plot_list[[4]] + plot_layout(ncol = 2,nrow = 2)
  setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
  ggsave(plot = cor_plots,filename=paste0('FDR1_',p_val,'_CorrPvalPlots.png'),width=6,height=4)
}

################################################################
#######  Investigate veQTL and their relation to eQTL ##########
################################################################

####################################################################################
# (1) How many significant gene-SNP veQTL pairs across the transcriptome-genome ####
####################################################################################

p_val_text <- paste0(' p < ',p_val)
print(paste0('Ctrl_cis_veQTL_sig_FDR1',p_val_text))
print(nrow(Ctrl_cis_veQTL_sig)) 
print(paste0('HS_cis_veQTL_sig_FDR1',p_val_text))
print(nrow(HS_cis_veQTL_sig)) 
print(paste0('Ctrl_trans_veQTL_sig_FDR1',p_val_text))
print(nrow(Ctrl_trans_veQTL_sig)) 
print(paste0('HS_trans_veQTL_sig_FDR1',p_val_text))
print(nrow(HS_trans_veQTL_sig)) 

print(paste0('Ctrl_cis_eqtl_sig_FDR1',p_val_text))
print(nrow(Ctrl_cis_eqtl_sig)) 
print(paste0('HS_cis_eqtl_sig_FDR1',p_val_text))
print(nrow(HS_cis_eqtl_sig)) 
print(paste0('Ctrl_trans_eqtl_sig_FDR1',p_val_text))
print(nrow(Ctrl_trans_eqtl_sig))
print(paste0('HS_trans_eqtl_sig_FDR1',p_val_text))
print(nrow(HS_trans_eqtl_sig))


#########################################################################################
# (1b) Dividing these in control-only shared, HS-only shared, and shared - pval #########
#########################################################################################

# Initialize the output dataframe
row_labels <- c("cis_FDR<0.01", "cis_FDR<0.1", "cis_FDR<0.05", 
                "trans_FDR<0.01", "trans_FDR<0.1", "trans_FDR<0.05",
                "cisandtrans_FDR<0.01", "cisandtrans_FDR<0.1", "cisandtrans_FDR<0.05")

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
FDR_cutoffs <- c(0.01, 0.1, 0.05)
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
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\GxE")
write.csv(veQTL_GxE_table,'FDR1_veQTL_GxE_table.csv')

# Turn it into a percentage
veQTL_GxE_table_percent <- veQTL_GxE_table
for (n_pval in 0:2){
  firstrow <- (n_pval*3)+1
  lastrow <- firstrow + 2
  veQTL_GxE_table_percent[firstrow:lastrow,2:4] <- veQTL_GxE_table[firstrow:lastrow,2:4]*100/(sum(veQTL_GxE_table[firstrow,2:4]))
  veQTL_GxE_table_percent[firstrow:lastrow,5:7] <- veQTL_GxE_table[firstrow:lastrow,5:7]*100/(sum(veQTL_GxE_table[firstrow,5:7]))
}
veQTL_GxE_table_percent[,2:7] <- round(veQTL_GxE_table_percent[,2:7],digits = 2)
View(veQTL_GxE_table_percent)
write.csv(veQTL_GxE_table_percent,'FDR1_veQTL_GxE_table_percent.csv')
# Make a plot of the percentages
veQTL_GxE_table_percent_separated <- veQTL_GxE_table_percent %>%
  separate(SharedCutoff, into = c("cis_or_trans", "shared_cutoff"), sep = "_")
# Shared is the opposite direction actually, turn the signs around
veQTL_GxE_table_percent_separated$shared_cutoff <- 
  factor(veQTL_GxE_table_percent_separated$shared_cutoff,
         levels = c("FDR<0.1", "FDR<0.05", "FDR<0.01")
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
  
  # Save plot
  ggsave(filename = paste0('FDR1_',col_name, "_veQTL_dotplot.svg"),
         plot = p, width = 6, height = 4, dpi = 300)
}


##############################################################
##### Save files and working environment so far ##############
##############################################################

# Save necessary columns of top veQTL
plotting_columns <- c("GENE", "CHROM","POS","SNP","REF","ALT","P","vfdr","COR")
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
write.csv(HS_cis_veQTL_sig[,plotting_columns],"FDR1_HS_cis_veQTL_sig.csv")
write.csv(Ctrl_cis_veQTL_sig[,plotting_columns],"FDR1_Ctrl_cis_veQTL_sig.csv")
write.csv(Ctrl_trans_veQTL_sig[,plotting_columns],"FDR1_Ctrl_trans_veQTL_sig.csv")
write.csv(HS_trans_veQTL_sig[,plotting_columns],"FDR1_HS_trans_veQTL_sig.csv")


setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
save.image(file='5_Investigating_eQTL_and_veQTL_FDR1.RData')
