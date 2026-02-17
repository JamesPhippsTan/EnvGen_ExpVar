# Rerunning the cis-mappings but using the non-quantile normalised phenotypes
# If there is an impact, then the disparity should be reflected in the number 
# Of cis-eQTL 

# Last updated: 22/10/2025

##############################
##### Packages and Setup #####
##############################

rm(list = ls())

library(data.table)
library(reshape)
library(tibble)
library(tidyr)

# Functions
# Quantile normalisation function to turn each gene's distribution into a standard normal distribution
quant_norm <- function (v) {
  w <- rep(0, length(v))
  ranking <- rank(v, ties.method = "min")
  n <- length(unique(ranking))
  w[1:length(v) %in% ranking] <- qnorm(seq(from = 1/(n + 1), 
                                           to = n/(n + 1), by = 1/(n + 1)))
  return(w[ranking])
}

# Apply quantile normalisation to a whole-transcriptome count matrix
quant_norm_count_matrix <- function (count_matrix) {
  count_matrix_normalised <- count_matrix
  for (genes in colnames(count_matrix_normalised)){
    count_matrix_normalised[,genes] <- quant_norm(count_matrix[,genes])
  }
  return(count_matrix_normalised)
}

######################
##### Datasets #######
######################

# GRAMMARed phenotypes for Ctrl
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
setwd("GRAMMARed_phenotypes_Ctrl")
GRAMMAR_phenotypes_Ctrl_list <- list.files(pattern = "\\.txt$")
GRAMMARed_phenotypes_Ctrl_table <- lapply(GRAMMAR_phenotypes_Ctrl_list, read.table)

# GRAMMARed phenotypes for HS
setwd("../GRAMMARed_phenotypes_HS/")
GRAMMAR_phenotypes_HS_list <- list.files(pattern = "\\.txt$")
GRAMMARed_phenotypes_HS_table <- lapply(GRAMMAR_phenotypes_HS_list, read.table)

# Gene body locations from when the data was generated (downloaded from the droseqtl.org)
# These are necessary for cis-mapping
setwd("..")
Gene_body_locations <-fread('droseQTLorg_genes_positions.csv',header = T,select = c(1,2,3,6))
# Substitute X with 23, as this is how the SNP chromosomes are encoded in the vcf files
Gene_body_locations$Chr <- gsub('X','23',Gene_body_locations$Chr)

######################################################################
##### (1) Create a single GRAMMARed phenotypes dataframe #############
######################################################################

# Merge the files
GRAMMARed_phenotypes_Ctrl <- merge_recurse(GRAMMARed_phenotypes_Ctrl_table, by = "id")
GRAMMARed_phenotypes_Ctrl <- column_to_rownames(GRAMMARed_phenotypes_Ctrl,var = 'id')
GRAMMARed_phenotypes_HS <- merge_recurse(GRAMMARed_phenotypes_HS_table, by = "id")
GRAMMARed_phenotypes_HS <- column_to_rownames(GRAMMARed_phenotypes_HS,var = 'id')

#######################################################################
##### (2) Get summary statistics of distributions before and after QN #
#######################################################################

# Before quantile normal transform
GRAMMARed_phenotypes_Ctrl_orig <- GRAMMARed_phenotypes_Ctrl
GRAMMARed_phenotypes_HS_orig <- GRAMMARed_phenotypes_HS

# After quantile normal transform
GRAMMARed_phenotypes_Ctrl_quant <- quant_norm_count_matrix(GRAMMARed_phenotypes_Ctrl)
GRAMMARed_phenotypes_HS_quant <- quant_norm_count_matrix(GRAMMARed_phenotypes_HS)

# Get the correlation before and after transformation
BeforeAfterQuanNormCorr <- data.frame(gene=character(),
                                      Ctrl_Corr=numeric(),
                                      HS_Corr=numeric(),
                                      KS_Before=numeric())

for (gene in colnames(GRAMMARed_phenotypes_Ctrl)){
GeneCors <- data.frame(gene=gene,
                       Ctrl_Corr= cor.test(GRAMMARed_phenotypes_Ctrl_orig[,gene],GRAMMARed_phenotypes_Ctrl_quant[,gene])$estimate,
                       HS_Corr= cor.test(GRAMMARed_phenotypes_HS_orig[,gene],GRAMMARed_phenotypes_HS_quant[,gene])$estimate,
                       KS_Before=ks.test(GRAMMARed_phenotypes_HS_orig[,gene],GRAMMARed_phenotypes_Ctrl_orig[,gene])$p)
BeforeAfterQuanNormCorr <- rbind(BeforeAfterQuanNormCorr,GeneCors)
}
summary(BeforeAfterQuanNormCorr$Ctrl_Corr) # 85%-99%; 98% mean and median
summary(BeforeAfterQuanNormCorr$HS_Corr) # 83%-99%; 96% mean and median
summary(BeforeAfterQuanNormCorr$KS_Before) # 0-0.99; 0.055 mean

# The distributions of most genes change quite little before and after the normalisation, given the high correlation coefficient
# Therefore, QN still could have an effect on some genes

# If QN induces the result that HS has many more hits than Ctrl
# Then this should be reflected at the cis-veQTL level
# Seeing that the number of cis-veQTL is already very much higher than 

#######################################################################
##### (3) Format into three distinct phenotype files per condition ####
#######################################################################

#####################################
##### (3a) Ctrl population ##########
#####################################

# Transpose and add gene names as a column
Phen_GRAMMARed_Ctrl <- data.frame(Gene=rownames(t(GRAMMARed_phenotypes_Ctrl_orig)),t(GRAMMARed_phenotypes_Ctrl_orig),check.names = F)

# Same rationale as before
Phen_GRAMMARed_Ctrl_positioned <- merge(Gene_body_locations,Phen_GRAMMARed_Ctrl,all.y = T,by = 'Gene')
Phen_GRAMMARed_Ctrl_positioned[,1:4] <- Phen_GRAMMARed_Ctrl_positioned[,c(4,2,3,1)]
colnames(Phen_GRAMMARed_Ctrl_positioned)[1:4] <- c('chr','start','end','phenotype_id')
Phen_GRAMMARed_Ctrl_positioned <- Phen_GRAMMARed_Ctrl_positioned[order(Phen_GRAMMARed_Ctrl_positioned$chr,Phen_GRAMMARed_Ctrl_positioned$start),]
colnames(Phen_GRAMMARed_Ctrl_positioned)[1] <- '#chr'

# For cis-eQTL and -veQTL
Phen_GRAMMARed_Ctrl_cis <- Phen_GRAMMARed_Ctrl_positioned[complete.cases(Phen_GRAMMARed_Ctrl_positioned),]

#####################################
##### (3b) HS population ############
#####################################

# Transpose and add gene names as a column
Phen_GRAMMARed_HS <- data.frame(Gene=rownames(t(GRAMMARed_phenotypes_HS_orig)),t(GRAMMARed_phenotypes_HS_orig),check.names = F)

# Same rationale as before
Phen_GRAMMARed_HS_positioned <- merge(Gene_body_locations,Phen_GRAMMARed_HS,all.y = T,by = 'Gene')
Phen_GRAMMARed_HS_positioned[,1:4] <- Phen_GRAMMARed_HS_positioned[,c(4,2,3,1)]
colnames(Phen_GRAMMARed_HS_positioned)[1:4] <- c('chr','start','end','phenotype_id')
Phen_GRAMMARed_HS_positioned <- Phen_GRAMMARed_HS_positioned[order(Phen_GRAMMARed_HS_positioned$chr,Phen_GRAMMARed_HS_positioned$start),]
colnames(Phen_GRAMMARed_HS_positioned)[1] <- '#chr'

# For cis-eQTL and -veQTL
Phen_GRAMMARed_HS_cis <- Phen_GRAMMARed_HS_positioned[complete.cases(Phen_GRAMMARed_HS_positioned),]

#############################################################
##### (3) Save the datasets into .bed files #################
#############################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\noQN")

write.table(Phen_GRAMMARed_HS_cis,"HS_phenos_cis_mapping_noQN.bed",sep='\t',quote=F,row.names=F)
write.table(Phen_GRAMMARed_Ctrl_cis,"Ctrl_phenos_cis_mapping_noQN.bed",sep='\t',quote=F,row.names=F)

# Rerun the cis-eQTL and cis-veQTL pipelines
# Assume the cis-eQTL used to map are the same 

######################################################################
##### (4) Check if cis-veQTL mapping results are very different ######
######################################################################

# Assuming that the cis-eQTL called are roughly the same 
# The eQTLs that go into the -eQTL tag should be the same

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")

add_SNP_col <- function(veQTL_df){
  veQTL_df$SNP <- paste(veQTL_df$CHROM,veQTL_df$POS)
return(veQTL_df)
}

# Cis-veQTL with QN
Ctrl_cis_veQTL <- add_SNP_col(read.table("Ctrl_cis_veQTL.txt", header = T,skipNul = T))
HS_cis_veQTL <- add_SNP_col(read.table("HS_cis_veQTL.txt", header = T,skipNul = T))
# Only difference is the lack of transformation

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\noQN/")

# Cis-veQTL without QN
Ctrl_cis_veQTL_noQN <- add_SNP_col(read.table("Ctrl_cis_veQTL_noQN.txt", header = T,skipNul = T))
HS_cis_veQTL_noQN <- add_SNP_col(read.table("HS_cis_veQTL_noQN.txt", header = T,skipNul = T))
# Only difference is the lack of transformation

n_cis_tests <- nrow(Ctrl_cis_veQTL)

# FDR correction

# With QN
Ctrl_cis_veQTL$vfdr <- p.adjust(Ctrl_cis_veQTL$P,method = 'BH',n=n_cis_tests)
Ctrl_cis_veQTL_sig <- subset(Ctrl_cis_veQTL,vfdr<0.05) 
max(Ctrl_cis_veQTL_sig$P) # 7.20328e-07
nrow(Ctrl_cis_veQTL_sig) # 11
length(unique(Ctrl_cis_veQTL_sig$GENE)) # 9 genes
length(unique(Ctrl_cis_veQTL_sig$SNP)) # 9 SNPs
HS_cis_veQTL$vfdr <- p.adjust(HS_cis_veQTL$P,method = 'BH',n=n_cis_tests)
HS_cis_veQTL_sig <- subset(HS_cis_veQTL,vfdr<0.05) 
max(HS_cis_veQTL_sig$P) # 1.00678e-05
nrow(HS_cis_veQTL_sig) # 148
length(unique(HS_cis_veQTL_sig$GENE)) # 130 genes
length(unique(HS_cis_veQTL_sig$SNP)) # 145 SNPs

# No QN
Ctrl_cis_veQTL_noQN$vfdr <- p.adjust(Ctrl_cis_veQTL_noQN$P,method = 'BH',n=n_cis_tests)
Ctrl_cis_veQTL_noQN_sig <- subset(Ctrl_cis_veQTL_noQN,vfdr<0.05) 
max(Ctrl_cis_veQTL_noQN_sig$P) # 1.60384e-05
nrow(Ctrl_cis_veQTL_noQN_sig) # 235
length(unique(Ctrl_cis_veQTL_noQN_sig$GENE)) # 114 genes
length(unique(Ctrl_cis_veQTL_noQN_sig$SNP)) # 221 SNPs
HS_cis_veQTL_noQN$vfdr <- p.adjust(HS_cis_veQTL_noQN$P,method = 'BH',n=n_cis_tests)
HS_cis_veQTL_noQN_sig <- subset(HS_cis_veQTL_noQN,vfdr<0.05) 
max(HS_cis_veQTL_noQN_sig$P) # 7.16092e-05
nrow(HS_cis_veQTL_noQN_sig) # 1049
length(unique(HS_cis_veQTL_noQN_sig$GENE)) # 586 genes
length(unique(HS_cis_veQTL_noQN_sig$SNP)) # 1006 SNPs
# No quantile normalisation dramatically increases the number of cis-veQTL using the BH method

# How many of the pairs intersect?
# Control
length(intersect(paste(Ctrl_cis_veQTL_sig$GENE,Ctrl_cis_veQTL_sig$SNP),
                 paste(Ctrl_cis_veQTL_noQN_sig$GENE,Ctrl_cis_veQTL_noQN_sig$SNP)))
# All 11 of the 11
# HS
length(intersect(paste(HS_cis_veQTL_sig$GENE,HS_cis_veQTL_sig$SNP),
                 paste(HS_cis_veQTL_noQN_sig$GENE,HS_cis_veQTL_noQN_sig$SNP)))
# 142 out of the 148 -> close but not perfect


# How many of the genes and SNPs intersect for HS?
length(intersect(unique(HS_cis_veQTL_sig$GENE),
                 unique(HS_cis_veQTL_noQN_sig$GENE)))
# 125
length(intersect(unique(HS_cis_veQTL_sig$SNP),
                 unique(HS_cis_veQTL_noQN_sig$SNP)))
# 139

# Takeaway
# Not running quantile normalisation 
# Still results in more HS than Ctrl veQTL
# The same overlapping set of top veQTL
# And simply increases the number of veQTL called  
# Therefore, QN is not responsible for the imbalance nor the high number of veQTL
# It excacerbates it
# Given the similar FDR approach and the increased number of tests
# We expect the same excacerbations to occur in the trans-veQTL mapping
# Resulting in an even more exagerated result if QN was not applied


# Try the GxE calling approach
row_labels <- c("cis_FDR<0.2", "cis_FDR<0.1", "cis_FDR<0.05")
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
FDR_cutoffs <- c(0.2, 0.1, 0.05)
for (FDR in FDR_cutoffs) {
  row_name <- paste0("cis_FDR<", FDR)
  
  Ctrl_cis_cutoff <- subset(Ctrl_cis_veQTL_noQN, vfdr < FDR)
  HS_cis_cutoff <- subset(HS_cis_veQTL_noQN, vfdr < FDR)
  
  Shared_genes <- length(unique(c(
    intersect(Ctrl_cis_veQTL_noQN_sig$GENE, HS_cis_cutoff$GENE),
    intersect(HS_cis_veQTL_noQN_sig$GENE, Ctrl_cis_cutoff$GENE)
  )))
  
  GxE_Ctrl_genes <- length(setdiff(Ctrl_cis_veQTL_noQN_sig$GENE, HS_cis_cutoff$GENE))
  GxE_HS_genes   <- length(setdiff(HS_cis_veQTL_noQN_sig$GENE, Ctrl_cis_cutoff$GENE))
  
  Shared_SNPs <- length(unique(c(
    intersect(Ctrl_cis_veQTL_noQN_sig$SNP, HS_cis_cutoff$SNP),
    intersect(HS_cis_veQTL_noQN_sig$SNP, Ctrl_cis_cutoff$SNP)
  )))
  
  GxE_Ctrl_SNPs <- length(setdiff(Ctrl_cis_veQTL_noQN_sig$SNP, HS_cis_cutoff$SNP))
  GxE_HS_SNPs   <- length(setdiff(HS_cis_veQTL_noQN_sig$SNP, Ctrl_cis_cutoff$SNP))
  
  veQTL_GxE_table[veQTL_GxE_table$SharedCutoff == row_name, 2:7] <- c(
    GxE_Ctrl_genes, GxE_HS_genes, Shared_genes,
    GxE_Ctrl_SNPs, GxE_HS_SNPs, Shared_SNPs
  )
}
View(veQTL_GxE_table)

# Plot as percentages
veQTL_GxE_table_percent <- veQTL_GxE_table
veQTL_GxE_table_percent[1:3,2:4] <- veQTL_GxE_table[1:3,2:4]*100/(sum(veQTL_GxE_table[1,2:4]))
veQTL_GxE_table_percent[1:3,5:7] <- veQTL_GxE_table[1:3,5:7]*100/(sum(veQTL_GxE_table[1,5:7]))
veQTL_GxE_table_percent[,2:7] <- round(veQTL_GxE_table_percent[,2:7],digits = 2)
View(veQTL_GxE_table_percent)


#####################################################################
##### (5) Check if cis-eQTL mapping results are very different ######
#####################################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")

# Cis-eQTL with QN
Ctrl_cis_eQTL <- read.table("Ctrl_cis_eqtl_result/Ctrl.cis_qtl.txt.gz", header = T,skipNul = T)
HS_cis_eQTL <- read.table("HS_cis_eqtl_result/HS.cis_qtl.txt.gz", header = T,skipNul = T)

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\noQN")

# Cis-eQTL without QN
Ctrl_cis_eQTL_noQN <- read.table("Ctrl_cis_eqtl_result_noQN/Ctrl_noQN.cis_qtl.txt.gz", header = T,skipNul = T)
HS_cis_eQTL_noQN <- read.table("HS_cis_eqtl_result_noQN/HSnoQN.cis_qtl.txt.gz", header = T,skipNul = T)
# Only difference is the lack of transformation

n_cis_tests <- nrow(Ctrl_cis_eQTL)

# FDR correction

# With QN
Ctrl_cis_eQTL$vfdr <- p.adjust(Ctrl_cis_eQTL$pval_nominal,method = 'BH',n=n_cis_tests)
Ctrl_cis_eQTL_sig <- subset(Ctrl_cis_eQTL,vfdr<0.05) 
nrow(Ctrl_cis_eQTL_sig) # 86250
length(unique(Ctrl_cis_eQTL_sig$phenotype_id)) # 4678 genes
length(unique(Ctrl_cis_eQTL_sig$variant_id)) # 71870 SNPs

HS_cis_eQTL$vfdr <- p.adjust(HS_cis_eQTL$pval_nominal,method = 'BH',n=n_cis_tests)
HS_cis_eQTL_sig <- subset(HS_cis_eQTL,vfdr<0.05) 
nrow(HS_cis_eQTL_sig) # 60189
length(unique(HS_cis_eQTL_sig$phenotype_id)) # 130 genes
length(unique(HS_cis_eQTL_sig$variant_id)) # 145 SNPs

# No QN
Ctrl_cis_eQTL_noQN$vfdr <- p.adjust(Ctrl_cis_eQTL_noQN$pval_nominal,method = 'BH',n=n_cis_tests)
Ctrl_cis_eQTL_noQN_sig <- subset(Ctrl_cis_eQTL_noQN,vfdr<0.05) 
nrow(Ctrl_cis_eQTL_noQN_sig) # 77445
length(unique(Ctrl_cis_eQTL_noQN_sig$phenotype_id)) # 4424 genes
length(unique(Ctrl_cis_eQTL_noQN_sig$variant_id)) # 65185 SNPs

HS_cis_eQTL_noQN$vfdr <- p.adjust(HS_cis_eQTL_noQN$pval_nominal,method = 'BH',n=n_cis_tests)
HS_cis_eQTL_noQN_sig <- subset(HS_cis_eQTL_noQN,vfdr<0.05) 
nrow(HS_cis_eQTL_noQN_sig) # 46968
length(unique(HS_cis_eQTL_noQN_sig$phenotype_id)) # 3137 genes
length(unique(HS_cis_eQTL_noQN_sig$variant_id)) # 1006 SNPs

# How many of the pairs intersect?
# Control
length(intersect(paste(Ctrl_cis_eQTL_sig$phenotype_id,Ctrl_cis_eQTL_sig$variant_id),
                 paste(Ctrl_cis_eQTL_noQN_sig$phenotype_id,Ctrl_cis_eQTL_noQN_sig$variant_id)))
# 75119
# HS
length(intersect(paste(HS_cis_eQTL_sig$phenotype_id,HS_cis_eQTL_sig$variant_id),
                 paste(HS_cis_eQTL_noQN_sig$phenotype_id,HS_cis_eQTL_noQN_sig$variant_id)))
# 45633


# How many of the genes and SNPs intersect for Ctrl?
length(intersect(unique(Ctrl_cis_eQTL_sig$phenotype_id),
                 unique(Ctrl_cis_eQTL_noQN_sig$phenotype_id)))
# 4277
length(intersect(unique(Ctrl_cis_eQTL_sig$variant_id),
                 unique(Ctrl_cis_eQTL_noQN_sig$variant_id)))
# 63430

length(intersect(unique(HS_cis_eQTL_sig$phenotype_id),
                 unique(HS_cis_eQTL_noQN_sig$phenotype_id)))
# 3017
length(intersect(unique(HS_cis_eQTL_sig$variant_id),
                 unique(HS_cis_eQTL_noQN_sig$variant_id)))
# 139

# Takeaway
# Not running quantile normalisation 
# Still results in fewer HS than Ctrl eQTL
# The same overlapping set of top eQTL
# And simply decreases the number of eQTL called  
# With fewer eQTL being controlled, this would increase the number of veQTL called further 

# Try the GxE calling approach
row_labels <- c("cis_FDR<0.2", "cis_FDR<0.1", "cis_FDR<0.05")
eQTL_GxE_table <- data.frame(
  SharedCutoff = row_labels,
  GxE_Ctrl_genes = integer(length(row_labels)),
  GxE_HS_genes = integer(length(row_labels)),
  Shared_genes = integer(length(row_labels)),
  GxE_Ctrl_SNPs = integer(length(row_labels)),
  GxE_HS_SNPs = integer(length(row_labels)),
  Shared_SNPs = integer(length(row_labels)),
  stringsAsFactors = FALSE
)
FDR_cutoffs <- c(0.2, 0.1, 0.05)
for (FDR in FDR_cutoffs) {
  row_name <- paste0("cis_FDR<", FDR)
  
  Ctrl_cis_cutoff <- subset(Ctrl_cis_eQTL_noQN, vfdr < FDR)
  HS_cis_cutoff <- subset(HS_cis_eQTL_noQN, vfdr < FDR)
  
  Shared_genes <- length(unique(c(
    intersect(Ctrl_cis_eQTL_noQN_sig$phenotype_id, HS_cis_cutoff$phenotype_id),
    intersect(HS_cis_eQTL_noQN_sig$phenotype_id, Ctrl_cis_cutoff$phenotype_id)
  )))
  
  GxE_Ctrl_genes <- length(setdiff(Ctrl_cis_eQTL_noQN_sig$phenotype_id, HS_cis_cutoff$phenotype_id))
  GxE_HS_genes   <- length(setdiff(HS_cis_eQTL_noQN_sig$phenotype_id, Ctrl_cis_cutoff$phenotype_id))
  
  Shared_SNPs <- length(unique(c(
    intersect(Ctrl_cis_eQTL_noQN_sig$variant_id, HS_cis_cutoff$variant_id),
    intersect(HS_cis_eQTL_noQN_sig$variant_id, Ctrl_cis_cutoff$variant_id)
  )))
  
  GxE_Ctrl_SNPs <- length(setdiff(Ctrl_cis_eQTL_noQN_sig$variant_id, HS_cis_cutoff$variant_id))
  GxE_HS_SNPs   <- length(setdiff(HS_cis_eQTL_noQN_sig$variant_id, Ctrl_cis_cutoff$variant_id))
  
  eQTL_GxE_table[eQTL_GxE_table$SharedCutoff == row_name, 2:7] <- c(
    GxE_Ctrl_genes, GxE_HS_genes, Shared_genes,
    GxE_Ctrl_SNPs, GxE_HS_SNPs, Shared_SNPs
  )
}

View(veQTL_GxE_table)
View(eQTL_GxE_table)


########################################


# Cis-veQTL without QN and without eQTL
Ctrl_cis_veQTL_noQN_noeQTL <- read.table("Ctrl_cis_veQTL_noQNnoeQTL.txt", header = T,skipNul = T)
HS_cis_veQTL_noQN_noeQTL <- read.table("HS_cis_veQTL_noQNnoeQTL.txt", header = T,skipNul = T)
# Only difference is the lack of transformation and the lack of eQTL regression

# Extra: No QN and No eQTL
Ctrl_cis_veQTL_noQN_noeQTL$vfdr <- p.adjust(Ctrl_cis_veQTL_noQN_noeQTL$pval_nominal,method = 'BH',n=n_cis_tests)
Ctrl_cis_veQTL_noQN_noeQTL_sig <- subset(Ctrl_cis_veQTL_noQN_noeQTL,vfdr<0.05) 
max(Ctrl_cis_veQTL_noQN_noeQTL_sig$pval_nominal) # 4.01684e-05
nrow(Ctrl_cis_veQTL_noQN_noeQTL_sig) # 589
length(unique(Ctrl_cis_veQTL_noQN_noeQTL_sig$phenotype_id)) # 244 genes
HS_cis_veQTL_noQN_noeQTL$vfdr <- p.adjust(HS_cis_veQTL_noQN_noeQTL$pval_nominal,method = 'BH',n=n_cis_tests)
HS_cis_veQTL_noQN_noeQTL_sig <- subset(HS_cis_veQTL_noQN_noeQTL,vfdr<0.05) 
max(HS_cis_veQTL_noQN_noeQTL_sig$pval_nominal) # 0.000114626
nrow(HS_cis_veQTL_noQN_noeQTL_sig) # 1678
length(unique(HS_cis_veQTL_noQN_noeQTL_sig$phenotype_id)) # 784 genes
# No eQTL tag dramatically increases the number of cis-veQTL using the BH method

# How many intersect? Both vs noQN and noeQTL
# Control
length(intersect(paste(Ctrl_cis_veQTL_sig$GENE,Ctrl_cis_veQTL_sig$POS),
                 paste(Ctrl_cis_veQTL_noQN_noeQTL_sig$GENE,Ctrl_cis_veQTL_noQN_noeQTL_sig$POS)))
# 9 of the 11
length(intersect(paste(HS_cis_veQTL_sig$GENE,HS_cis_veQTL_sig$POS),
                 paste(HS_cis_veQTL_noQN_noeQTL_sig$GENE,HS_cis_veQTL_noQN_noeQTL_sig$POS)))
# 144 out of the 148
# 90% of the QN eQTL-controlled pairs are within the
# significant pairs of mappings on QN eQTL-controlled or QN eQTL non-controlled data
# Which means the transformations are indeed more stringent and not looser

