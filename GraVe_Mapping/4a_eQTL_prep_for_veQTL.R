# Preparing eQTL to be regressed out of model to detect veQTL
# The idea is that some veQTL are caused by interacting eQTL variants having LD with a neighbouring eQTL
# The way to address this is to regress the effect of the interacting eQTL variant
# Only 50 can be mapped at maximum on our cluster...any more and issues arise

# Note: big bug in the top_eQTL_list function discovered 
# Set as decreasing = TRUE for the initial run 
# The correct parameter is decreasing = FALSE
# This has been corrected and the code has been rerun


# Last Updated: 27/3/26

#################################
##### Packages and Setup ########
#################################

rm(list = ls())
library(tibble)
library(stats)
library(data.table)
library(UpSetR)
library(ggplot2)
library(ggvenn)

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
load(file='4a_eQTL_prep_for_veQTL.RData')

# Functions
eQTL_number_per_gene <- function(gene_df){
  split_gene_df <- split(gene_df,f=gene_df[,'phenotype_id'])
  number_of_eQTL_by_gene <- as.data.frame(do.call(rbind,lapply(split_gene_df, nrow)))
  number_of_eQTL_by_gene$phenotype_id <- rownames(number_of_eQTL_by_gene)
  return(number_of_eQTL_by_gene)
}

top_eQTL_list <- function(split_gene_df,number) {
  sorted_df <- split_gene_df[order(split_gene_df[, 'p_final'], decreasing = FALSE), ]
  top_list <- head(sorted_df, number)
  return(top_list)
}

eQTL_correct_format <- function(eqtl_list) {
  formatted_eqtl_list <- eqtl_list[,c(7,2,3,5,6)] # corresponds to gene SNP_chrom SNP_pos allele0 allele1
  colnames(formatted_eqtl_list) <- c("gene", "SNP_chrom", "SNP_pos", "allele0", "allele1")
  return(formatted_eqtl_list)
}

######################
##### Datasets #######
######################

# Obtain eQTL results - results of all cis tests and all trans tests with pvals < 0.005
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
Ctrl_cis_eqtl_all <- read.table('Ctrl_cis_eqtl_result/Ctrl.cis_qtl.txt.gz',header = T, fill=TRUE)
HS_cis_eqtl_all <- read.table('HS_cis_eqtl_result/HS.cis_qtl.txt.gz',header = T, fill=TRUE)
Ctrl_trans_eqtl_all <- read.table('Ctrl_trans_eqtl_result/Ctrl.trans_qtl_pairs.txt.gz',header = T, fill=TRUE)
HS_trans_eqtl_all <- read.table('HS_trans_eqtl_result/HS.trans_qtl_pairs.txt.gz',header = T, fill=TRUE)

# Obtain eQTL positions
SNPs_dummy_positions <- read.table('SNPs_dummy_positions.txt')
SNPs_original_positions <-read.table('SNPs_original_positions.txt')


######################################################################
##### (1) Call significant eQTL using FDR correction procedure #######
######################################################################

# For cis tests (old): Take the beta p-values, multiply by number of genes, 0.05 cutoff
# For cis tests (new) and trans tests: Take the raw p-values, apply BH procedure to control FDR at 0.05

# Cis tests - BH method
n_cis_tests <- nrow(Ctrl_cis_eqtl_all)
Ctrl_cis_eqtl_all$p_final <- p.adjust(Ctrl_cis_eqtl_all$pval_nominal,method = 'BH',n = n_cis_tests)
HS_cis_eqtl_all$p_final <- p.adjust(HS_cis_eqtl_all$pval_nominal,method = 'BH',n = n_cis_tests)

# Cis tests - Permutation method combined with Bonferroni test
n_cis_tests <- nrow(Ctrl_cis_eqtl_all)
Ctrl_cis_eqtl_all$p_final_perm <- Ctrl_cis_eqtl_all$pval_beta_max*8763
HS_cis_eqtl_all$p_final_perm <- HS_cis_eqtl_all$pval_beta_max*8763

# Trans tests - BH method
n_trans_tests <- 8763*383710 - n_cis_tests
Ctrl_trans_eqtl_all$p_final <- p.adjust(Ctrl_trans_eqtl_all$pval,method = 'BH',n = n_trans_tests)
HS_trans_eqtl_all$p_final <- p.adjust(HS_trans_eqtl_all$pval,method = 'BH',n = n_trans_tests)

# Adjusted p-values threshold
p_val = 0.05

# Control and HS significant SNP set -> cis
Ctrl_cis_eqtl_sig <- subset(Ctrl_cis_eqtl_all,p_final<p_val)
HS_cis_eqtl_sig <- subset(HS_cis_eqtl_all,p_final<p_val)
nrow(Ctrl_cis_eqtl_sig)
nrow(HS_cis_eqtl_sig)

# Permutation verison
Ctrl_cis_eqtl_sig_perm <- subset(Ctrl_cis_eqtl_all,p_final_perm<p_val)
HS_cis_eqtl_sig_perm <- subset(HS_cis_eqtl_all,p_final_perm<p_val)

# Control and HS significant SNP set -> trans
Ctrl_trans_eqtl_sig <- subset(Ctrl_trans_eqtl_all,p_final<p_val)
max(Ctrl_trans_eqtl_sig$pval)
HS_trans_eqtl_sig <- subset(HS_trans_eqtl_all,p_final<p_val)
max(HS_trans_eqtl_sig$pval)
  
# Absolute effect size vs pval distributions as a check
plot(abs(Ctrl_cis_eqtl_sig$slope),Ctrl_cis_eqtl_sig$pval_beta_max)
plot(abs(HS_cis_eqtl_sig$slope),HS_cis_eqtl_sig$pval_beta_max)
plot(abs(Ctrl_trans_eqtl_sig$b),Ctrl_trans_eqtl_sig$pval)
plot(abs(HS_trans_eqtl_sig$b),HS_trans_eqtl_sig$pval)

# How many significant eQTL per gene
Ctrl_cis_eqtl_sig_number <- eQTL_number_per_gene(Ctrl_cis_eqtl_sig) 
HS_cis_eqtl_sig_number <- eQTL_number_per_gene(HS_cis_eqtl_sig) 

Ctrl_cis_eqtl_sig_number_perm <- eQTL_number_per_gene(Ctrl_cis_eqtl_sig_perm) 
HS_cis_eqtl_sig_number_perm <- eQTL_number_per_gene(HS_cis_eqtl_sig_perm) 

Ctrl_trans_eqtl_sig_number <- eQTL_number_per_gene(Ctrl_trans_eqtl_sig) 
HS_trans_eqtl_sig_number <- eQTL_number_per_gene(HS_trans_eqtl_sig) 


#######################################################
##### (2) Retain a fixed number of SNPs per gene ######
#######################################################

# Get columns needed
Cols_needed <- c('phenotype_id','variant_id','p_final')

############################
##### (2a) cis-veQTL #######
############################

# For cis-veQTL, eQTL in LD will likely be located within the cis-window 
# Use perm values for 'extra' certainty that the SNPs are causal
# Split into SNP by gene
Ctrl_cis_eQTL_by_gene <- split(Ctrl_cis_eqtl_sig_perm[,Cols_needed],f=Ctrl_cis_eqtl_sig_perm[,'phenotype_id'])
HS_cis_eQTL_by_gene <- split(HS_cis_eqtl_sig_perm[,Cols_needed],f=HS_cis_eqtl_sig_perm[,'phenotype_id'])

# Retain only the top 50 lowest-pval cis-eQTL 
Ctrl_top_cis_eQTL_list <- do.call(rbind,lapply(Ctrl_cis_eQTL_by_gene, top_eQTL_list,number=50))
HS_top_cis_eQTL_list <- do.call(rbind,lapply(HS_cis_eQTL_by_gene, top_eQTL_list,number=50))
# See below

# Number of SNPs per gene - just to get an idea of the distribution of high-eQTL genes
Ctrl_cis_eQTL_by_gene_number <- do.call(rbind,lapply(Ctrl_cis_eQTL_by_gene, nrow))
View(Ctrl_cis_eQTL_by_gene_number) 
# 29 genes below 50 cis-eQTL - 0.3% of all genes with expression
HS_cis_eQTL_by_gene_number <- do.call(rbind,lapply(HS_cis_eQTL_by_gene, nrow))
View(HS_cis_eQTL_by_gene_number) 
# 17 genes below 50 cis-eQTL - 0.2% of all genes with expression
# Imperfect eQTL accounting should not affect most of the transcriptome

############################
##### (2b) trans-veQTL #####
############################

# For trans-veQTL, the LD could be with variants anywhere, so all eQTL should be considered
# Concatenate cis and trans per condition
Ctrl_eqtl_sig <- rbind(Ctrl_cis_eqtl_sig[,Cols_needed],Ctrl_trans_eqtl_sig[,Cols_needed])
HS_eqtl_sig <- rbind(HS_cis_eqtl_sig[,Cols_needed],HS_trans_eqtl_sig[,Cols_needed])
nrow(Ctrl_eqtl_sig) 
nrow(HS_eqtl_sig) 

# Split into SNP by gene
Ctrl_eQTL_by_gene <- split(Ctrl_eqtl_sig,f=Ctrl_eqtl_sig[,'phenotype_id'])
HS_eQTL_by_gene <- split(HS_eqtl_sig,f=HS_eqtl_sig[,'phenotype_id'])

# Retain only the top 50 lowest-pval eQTL 
Ctrl_top_eQTL_list <- do.call(rbind,lapply(Ctrl_eQTL_by_gene, top_eQTL_list,number=50))
HS_top_eQTL_list <- do.call(rbind,lapply(HS_eQTL_by_gene, top_eQTL_list,number=50))
# See below note

# Number of SNPs per gene - just to get an idea of the distribution of high-eQTL genes
Ctrl_eQTL_by_gene_number <- do.call(rbind,lapply(Ctrl_eQTL_by_gene, nrow))
View(Ctrl_eQTL_by_gene_number)  
# 1434 genes above 50 eQTL - 16% of all genes with expression
HS_eQTL_by_gene_number <- do.call(rbind,lapply(HS_eQTL_by_gene, nrow))
View(HS_eQTL_by_gene_number)  
# 909 genes above 50 eQTL - 10% of all genes with expression
# Imperfect eQTL accounting should not affect most of the transcriptome

###########################################################
##### (4) Format list of top eQTL for veQTL mapping #######
###########################################################

# Rename the position column
SNPs_dummy_positions$variant_id <- SNPs_dummy_positions$V3
SNPs_original_positions$variant_id <- SNPs_original_positions$V3

# Find the positions of significant SNPs
# for cis
Ctrl_eqtl_sig_cis <- merge(SNPs_original_positions,Ctrl_top_cis_eQTL_list,all.y=T,by = 'variant_id')
HS_eqtl_sig_cis <- merge(SNPs_original_positions,HS_top_cis_eQTL_list,all.y=T,by = 'variant_id')
# for trans
Ctrl_eqtl_sig_trans <- merge(SNPs_dummy_positions,Ctrl_top_eQTL_list,all.y=T,by = 'variant_id')
HS_eqtl_sig_trans <- merge(SNPs_dummy_positions,HS_top_eQTL_list,all.y=T,by = 'variant_id')

######################################################
##### End) Save properly-formatted head eQTLs ########
######################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")

# Full significant eQTL list by condition and cis/trans
write.table(Ctrl_cis_eqtl_sig,'Ctrl_cis_eqtl_sig.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(HS_cis_eqtl_sig,'HS_cis_eqtl_sig.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(Ctrl_cis_eqtl_sig_perm,'Ctrl_cis_eqtl_sig_perm.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(HS_cis_eqtl_sig_perm,'HS_cis_eqtl_sig_perm.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(Ctrl_trans_eqtl_sig,'Ctrl_trans_eqtl_sig.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(HS_trans_eqtl_sig,'HS_trans_eqtl_sig.txt',row.names = F,col.names = T, quote = F, sep = " ")

# Number of eQTL per gene by condition and cis/trans
write.table(Ctrl_cis_eqtl_sig_number,'Ctrl_cis_eqtl_sig_number.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(HS_cis_eqtl_sig_number,'HS_cis_eqtl_sig_number.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(Ctrl_cis_eqtl_sig_number_perm,'Ctrl_cis_eqtl_sig_number_perm.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(HS_cis_eqtl_sig_number_perm,'HS_cis_eqtl_sig_number_perm.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(Ctrl_trans_eqtl_sig_number,'Ctrl_trans_eqtl_sig_number.txt',row.names = F,col.names = T, quote = F, sep = " ")
write.table(HS_trans_eqtl_sig_number,'HS_trans_eqtl_sig_number.txt',row.names = F,col.names = T, quote = F, sep = " ")

# Top 50 eQTL to be regressed out during veQTL mapper
# Old permutation method
write.table(eQTL_correct_format(Ctrl_eqtl_sig_cis),'eQTL_for_Ctrl_cis_veqtl_correct.txt',row.names = F,col.names = F, quote = F, sep = " ")
write.table(eQTL_correct_format(HS_eqtl_sig_cis),'eQTL_for_HS_cis_veqtl_correct.txt',row.names = F,col.names = F, quote = F, sep = " ")
write.table(eQTL_correct_format(Ctrl_eqtl_sig_trans),'eQTL_for_Ctrl_trans_veqtl_correct.txt',row.names = F,col.names = F, quote = F, sep = " ")
write.table(eQTL_correct_format(HS_eqtl_sig_trans),'eQTL_for_HS_trans_veqtl_correct.txt',row.names = F,col.names = F, quote = F, sep = " ")

save.image(file='4_eQTL_prep_for_veQTL.RData')

##################################################
##### Note on a bug and the extent of damage #####
################################################## 

# Triple checking if the top_eQTL_list() function works as intended - pick a gene with >50, see if it is correct
View(Ctrl_cis_eQTL_by_gene)
View(Ctrl_cis_eQTL_by_gene["FBgn0000479"][["FBgn0000479"]])
View(top_eQTL_list(Ctrl_cis_eQTL_by_gene["FBgn0000479"][["FBgn0000479"]],50))
View(subset(Ctrl_top_cis_eQTL_list,phenotype_id=='FBgn0000479'))
# It is wrong, it is printing the 50 HIGHEST pvals, Huiting was right...
# Corrected in 23.3.26
View(Ctrl_eQTL_by_gene)
View(Ctrl_eQTL_by_gene["FBgn0030313"][["FBgn0030313"]]) # full list
View(top_eQTL_list(Ctrl_eQTL_by_gene["FBgn0030313"][["FBgn0030313"]],50)) # function
View(subset(Ctrl_top_eQTL_list,phenotype_id=='FBgn0030313')) # output list
# It is wrong, it is printing the 50 HIGHEST pvals, Huiting was right...

# How many genes are affected by this?
# Number of SNPs per gene - just to get an idea of the distribution of high-eQTL genes
# Trans-veQTL
Ctrl_eQTL_by_gene_number <- do.call(rbind,lapply(Ctrl_eQTL_by_gene, nrow))
View(Ctrl_eQTL_by_gene_number)  
# 1443 genes above 50 eQTL - 16% of all genes with expression
HS_eQTL_by_gene_number <- do.call(rbind,lapply(HS_eQTL_by_gene, nrow))
View(HS_eQTL_by_gene_number)  
# 908 genes above 50 eQTL - 10% of all genes with expression
# Imperfect eQTL accounting should not affect most of the transcriptome

# Rerun the cis-veQTL mapping with the same datafiles
# Rerun the trans-veQTL mapping with only the above 10% - 15% of genes
# Just for the trans results because this is time-consuming
Ctrl_eQTL_by_gene_number_df <- as.data.frame(Ctrl_eQTL_by_gene_number)
HS_eQTL_by_gene_number_df <- as.data.frame(HS_eQTL_by_gene_number)
Ctrl_trans_veQTL_genes_remap <- rownames(subset(Ctrl_eQTL_by_gene_number_df,V1>50))
HS_trans_veQTL_genes_remap <- rownames(subset(HS_eQTL_by_gene_number_df,V1>50))
write.table(Ctrl_trans_veQTL_genes_remap,"Ctrl_trans_veQTL_genes_to_remap.txt",sep='\t',quote=F,row.names=F)
write.table(HS_trans_veQTL_genes_remap,"HS_trans_veQTL_genes_to_remap.txt",sep='\t',quote=F,row.names=F)


# Reduce the Phenotype dataframes to just those genes
# Create hashless copies by removing the # in front of chr so that the column names 'IDs' can be read
Ctrl_trans_veQTL_Phen <- read.table("Ctrl_phenos_trans_veQTL_mapping_hashless.bed",sep='\t',header = T,)
Ctrl_trans_veQTL_Phen_remap <- Ctrl_trans_veQTL_Phen[Ctrl_trans_veQTL_Phen$phenotype_id %in% Ctrl_trans_veQTL_genes_remap, ]
colnames(Ctrl_trans_veQTL_Phen_remap)[1] <-'#chr'
colnames(Ctrl_trans_veQTL_Phen_remap) <- gsub('X','',colnames(Ctrl_trans_veQTL_Phen_remap))
write.table(Ctrl_trans_veQTL_Phen_remap,"Ctrl_phenos_remap_trans_veQTL_mapping.bed",sep='\t',quote=F,row.names=F)

HS_trans_veQTL_Phen <- read.table("HS_phenos_trans_veQTL_mapping_hashless.bed",sep='\t',header = T,)
HS_trans_veQTL_Phen_remap <- HS_trans_veQTL_Phen[HS_trans_veQTL_Phen$phenotype_id %in% HS_trans_veQTL_genes_remap, ]
colnames(HS_trans_veQTL_Phen_remap)[1] <-'#chr'
colnames(HS_trans_veQTL_Phen_remap) <- gsub('X','',colnames(HS_trans_veQTL_Phen_remap))
write.table(HS_trans_veQTL_Phen_remap,"HS_phenos_remap_trans_veQTL_mapping.bed",sep='\t',quote=F,row.names=F)
