# Investigating eQTL and veQTL

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
load(file='5_Investigating_eQTL_and_veQTL.RData')

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

# Effect size (correlation slope) distributions
cor_plot_list <- list()

cor_plot_list[[1]] <-ggplot(Ctrl_cis_veQTL_sig,aes(y=abs(COR),x=P))+geom_point()+theme_classic()+geom_vline(xintercept = Ctrl_cis_veQTL_FDR0.05, color = "red", linetype = "dashed")+
  geom_vline(xintercept = Ctrl_cis_veQTL_FDR0.1, color = "red", linetype = "dashed")+
  ggtitle('Ctrl_cis_veQTL_sig')
cor_plot_list[[2]] <-ggplot(HS_cis_veQTL_sig,aes(y=abs(COR),x=P))+geom_point() +theme_classic()+geom_vline(xintercept = HS_cis_veQTL_FDR0.05, color = "red", linetype = "dashed")+
  geom_vline(xintercept = HS_cis_veQTL_FDR0.1, color = "red", linetype = "dashed")+
  ggtitle('HS_cis_veQTL_sig')
cor_plot_list[[3]] <-ggplot(Ctrl_trans_veQTL_sig,aes(y=abs(COR),x=P))+geom_point() +theme_classic()+geom_vline(xintercept = Ctrl_trans_veQTL_FDR0.05, color = "red", linetype = "dashed")+
  geom_vline(xintercept = Ctrl_trans_veQTL_FDR0.1, color = "red", linetype = "dashed")+
  ggtitle('Ctrl_trans_veQTL_sig')
cor_plot_list[[4]] <-ggplot(HS_trans_veQTL_sig,aes(y=abs(COR),x=P))+geom_point() +theme_classic()+geom_vline(xintercept = HS_trans_veQTL_FDR0.05, color = "red", linetype = "dashed")+
  geom_vline(xintercept = HS_trans_veQTL_FDR0.1, color = "red", linetype = "dashed")+
  ggtitle('HS_trans_veQTL_sig')

cor_plots <- cor_plot_list[[1]] + cor_plot_list[[2]] + cor_plot_list[[3]] + cor_plot_list[[4]] + plot_layout(ncol = 2,nrow = 2)
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
ggsave(plot = cor_plots,filename=paste0(p_val,'_CorrPvalPlots.png'),width=6,height=4)
}

################################################################
#######  Investigate veQTL and their relation to eQTL ##########
################################################################

####################################################################################
# (1) How many significant gene-SNP veQTL pairs across the transcriptome-genome ####
####################################################################################

p_val_text <- paste0(' p < ',p_val)
print(paste0('Ctrl_cis_veQTL_sig',p_val_text))
print(nrow(Ctrl_cis_veQTL_sig)) # 11
print(paste0('HS_cis_veQTL_sig',p_val_text))
print(nrow(HS_cis_veQTL_sig)) # 148
print(paste0('Ctrl_trans_veQTL_sig',p_val_text))
print(nrow(Ctrl_trans_veQTL_sig)) # 724
print(paste0('HS_trans_veQTL_sig',p_val_text))
print(nrow(HS_trans_veQTL_sig)) # 531266

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
row_labels <- c("cis_FDR<0.2", "cis_FDR<0.1", "cis_FDR<0.05", 
                "trans_FDR<0.2", "trans_FDR<0.1", "trans_FDR<0.05",
                "cisandtrans_FDR<0.2", "cisandtrans_FDR<0.1", "cisandtrans_FDR<0.05")

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
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\GxE")
write.csv(veQTL_GxE_table,'veQTL_GxE_table.csv')

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
write.csv(veQTL_GxE_table_percent,'veQTL_GxE_table_percent.csv')
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
  
  # Save plot
  ggsave(filename = paste0(col_name, "_veQTL_dotplot.svg"),
         plot = p, width = 6, height = 4, dpi = 300)
}

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
  Value = c(8763-0-1246-6921, 0, 1246, 6921) 
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
ggsave(plot = eQTL_plot, filename='eGene_overlap_plot.svg', width = 3, height = 1.5, dpi = 300)
ggsave(plot = veQTL_plot, filename='vGene_overlap_plot.svg', width = 3, height =1.5, dpi = 300)

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
Genes_with_no_QTL <- setdiff(All_genes$gene_id,QTL_number_per_gene$gene_id)
# Create a new tibble with 0s for other columns
new_rows <- QTL_number_per_gene %>%
  tibble(gene_id = Genes_with_no_QTL) %>%
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
    ggsave(plot = mean_plots,filename = paste0("MeanCor_",QTL_category,"plot.svg"),width=3,height=6,dpi=300)
    var_plots <- wrap_plots(var_plot_list, nrow = 2)  # Arrange in 2 columns
    print(var_plots)
    ggsave(plot = var_plots,filename = paste0("VarCor_",QTL_category,"plot.svg"),width=3,height=6,dpi=300)
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

Non_veQTL_SNPs = 383.710-GxE_Ctrl_SNPs/1000-Shared_SNPs/1000-GxE_HS_SNPs/1000
data_veQTL_SNPs <- data.frame(
  Partition = c("Non-veQTL", "Ctrl-only veQTL", "Ctrl and HS veQTL", "HS-only veQTL"),
  Value = c(383.710-0.079-108.446-2.622, 0.079, 2.622,108.446) 
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
ggsave(plot = eQTL_SNP_plot, filename='eQTL_SNP_overlap_plot.svg', width = 3, height = 1.5, dpi = 300)
ggsave(plot = veQTL_SNP_plot, filename='veQTL_SNP_overlap_plot.svg', width = 3, height =1.5, dpi = 300)


##############################################################################
# (3b) How many genes per SNP? Are there hotspot SNPs that regulate many genes?
##############################################################################

# Get the number of hotspot QTL for each test within one table
hotspot_SNPs <- merge(data.frame(table(Ctrl_cis_eqtl_sig$variant_id)),data.frame(table(Ctrl_trans_eqtl_sig$variant_id)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(HS_cis_eqtl_sig$variant_id)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(HS_trans_eqtl_sig$variant_id)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(Ctrl_cis_veQTL_sig$SNP)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(HS_cis_veQTL_sig$SNP)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(Ctrl_trans_veQTL_sig$SNP)),all = T,by='Var1')
hotspot_SNPs <- merge(hotspot_SNPs,data.frame(table(HS_trans_veQTL_sig$SNP)),all = T,by='Var1')
colnames(hotspot_SNPs) <- c('SNP','Ctrl_cis_eQTLs','Ctrl_trans_eQTLs','HS_cis_eQTLs','HS_trans_eQTLs','Ctrl_cis_veQTLs','HS_cis_veQTLs','Ctrl_trans_veQTLs','HS_trans_veQTLs')
hotspot_SNPs_summary <- as.data.frame(t(sapply(hotspot_SNPs[,2:ncol(hotspot_SNPs)], summary_with_count)))
View(hotspot_SNPs_summary)

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
  ggsave(plot = eQTL_pleiotropy,filename = paste0("SNP_pleiotropy,",SNP_column,"_plot.svg"),height=3,width=3,dpi=300)
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
  ggsave(plot = veQTL_pleiotropy,filename = paste0("SNP_pleiotropy,",SNP_column,"_plot.svg"),height=3,width=3,dpi=300)
}

###############################################################################
# (3c) Do the SNP sets have biases in minor allele frequency? eQTL vs veQTL ###
###############################################################################

# MAFs 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\MAF")
# Contrl and HS overlapping density plots
MAF_long <- data.frame(MAF = c(MAF$Ctrl_MAF,MAF$HS_MAF),
                       Condition = c(rep("Ctrl", length(MAF$Ctrl_MAF)),
                         rep("HS", length(MAF$HS_MAF))))

# MAF density plot
MAF_density <- ggplot(MAF_long, aes(x = MAF, col = Condition)) +
  geom_density(adjust = 2) +
  scale_color_manual(values = c("Ctrl" = "#C6B49E", "HS" = "#DF9F65")) +  
  theme_classic() +
  labs(x = "MAF", y = "Density", fill = "Condition") +
  theme(legend.position = c(0.25, 0.25),            # Position legend inside plot
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))+coord_cartesian(ylim=c(0.9,4),xlim=c(0,0.5))
MAF_density
ggsave(plot = MAF_density,filename = paste0("HS_Ctrl_MAF_Density_plot.svg"),height=3,width=3,dpi=300)

# Control vs HS MAF
MAF_spearman_test <- cor.test(y=MAF$HS_MAF, x=MAF$Ctrl_MAF, method = "spearman")
MAF_spearman_coef <- MAF_spearman_test$estimate
MAF_spearman_p_value <- MAF_spearman_test$p.value
MAF_spearman_test_result <-paste0("rho = ",format_statistic(MAF_spearman_coef),
                                   "\n","p = ",format_P_handle0(MAF_spearman_p_value))

MAF_correlation <- ggplot(MAF, aes(y = HS_MAF, x = Ctrl_MAF)) +
  geom_point(size = 0.8, alpha = 0.5,col='darkgrey') +
  geom_smooth(method = "lm", col = "blue",linetype='dashed')+
  annotate("text", y = max(MAF$HS_MAF), x = min(MAF$HS_MAF), label = MAF_spearman_test_result, size = 3, hjust = 0,vjust = 1)+
  ylab('MAF (HS)')+
  xlab('MAF (Ctrl)')+
  theme_classic()+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
  plot.background = element_rect(fill = "transparent", color = NA))+
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 0.5)) +
  scale_x_continuous(breaks = c(0, seq(0.1, 0.5, 0.1))) +
  scale_y_continuous(breaks = c(0, seq(0.1, 0.5, 0.1)))
MAF_correlation
ggsave(plot = MAF_correlation,filename = paste0("HS_Ctrl_MAF_Correlation_plot.svg"),height=2,width=2, bg = "transparent",dpi=300)
# Note: This is a HUGE svg because it contains 400K points

# eQTL list
eQTL_list <- list(
  All = unique(c(Ctrl_cis_eqtl_sig$variant_id,
                      Ctrl_trans_eqtl_sig$variant_id,
                      HS_cis_eqtl_sig$variant_id,
                    HS_trans_eqtl_sig$variant_id)),
  All_cis = unique(c(Ctrl_cis_eqtl_sig$variant_id,
                         HS_cis_eqtl_sig$variant_id)),
  All_trans = unique(c(Ctrl_trans_eqtl_sig$variant_id,
                            HS_trans_eqtl_sig$variant_id))
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

compare_MAFs_eQTL_veQTL(eQTL_datasets=eQTL_list, 
                                           veQTL_datasets=veQTL_list, 
                                           MAF_df=MAF)

##############################################################
##### Save files and working environment so far ##############
##############################################################

# Save necessary columns of top veQTL
plotting_columns <- c("GENE", "CHROM","POS","SNP","REF","ALT","P","vfdr","COR")
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
write.csv(HS_cis_veQTL_sig[,plotting_columns],"HS_cis_veQTL_sig.csv")
write.csv(Ctrl_cis_veQTL_sig[,plotting_columns],"Ctrl_cis_veQTL_sig.csv")
write.csv(Ctrl_trans_veQTL_sig[,plotting_columns],"Ctrl_trans_veQTL_sig.csv")
write.csv(HS_trans_veQTL_sig[,plotting_columns],"HS_trans_veQTL_sig.csv")

# Save hotspot genes 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_Genes/")
# Number of eQTL and veQTL per gene (with at least 1 eQTL or veQTL)
write.csv(QTL_number_per_gene_final,'QTL_number_per_gene.csv')
# Number of QTL for eGenes and vGenes
write.csv(QTL_number_per_gene_summary,'QTL_number_per_gene_summary.csv')

# Save hotspot SNPs
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs/")
write.csv(hotspot_SNPs,'hotspot_SNPs.csv')
write.csv(hotspot_SNPs_summary,'hotspot_SNPs_summaries.csv')

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
save.image(file='5_Investigating_eQTL_and_veQTL.RData')

