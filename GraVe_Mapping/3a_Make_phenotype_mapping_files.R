# Making the .bed phenotype files from the GRAMMAR-corrected counts
# Works for both eQTL and veQTL mapping
# Both algorithms use the same input phenotype file formats but with different specifications

# Last updated: 11/12/2024

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

# Quantile normal transform
GRAMMARed_phenotypes_Ctrl <- quant_norm_count_matrix(GRAMMARed_phenotypes_Ctrl)
GRAMMARed_phenotypes_HS <- quant_norm_count_matrix(GRAMMARed_phenotypes_HS)

# Transpose and add gene names as a column
Phen_GRAMMARed_Ctrl <- data.frame(Gene=rownames(t(GRAMMARed_phenotypes_Ctrl)),t(GRAMMARed_phenotypes_Ctrl),check.names = F)
Phen_GRAMMARed_HS <- data.frame(Gene=rownames(t(GRAMMARed_phenotypes_HS)),t(GRAMMARed_phenotypes_HS),check.names = F)

#############################################################
##### (2) Format into three distinct phenotype files per condition
#############################################################

#####################################
##### (2a) Control population #######
#####################################

# For all, we need to add gene locations as the first four columns in the correct order
Phen_GRAMMARed_Ctrl_positioned <- merge(Gene_body_locations,Phen_GRAMMARed_Ctrl,all.y = T,by = 'Gene')
Phen_GRAMMARed_Ctrl_positioned[,1:4] <- Phen_GRAMMARed_Ctrl_positioned[,c(4,2,3,1)]
colnames(Phen_GRAMMARed_Ctrl_positioned)[1:4] <- c('chr','start','end','phenotype_id')

# Then we need to rearrange in ascending order and rename the first column to include the #
Phen_GRAMMARed_Ctrl_positioned <- Phen_GRAMMARed_Ctrl_positioned[order(Phen_GRAMMARed_Ctrl_positioned$chr,Phen_GRAMMARed_Ctrl_positioned$start),]
colnames(Phen_GRAMMARed_Ctrl_positioned)[1] <- '#chr'

# For cis-eQTL and -veQTL
# No mitochondrial SNPs are in the dataset therefore we can exclude mitochondrial genes 
# One file works for both
Phen_GRAMMARed_Ctrl_cis <- Phen_GRAMMARed_Ctrl_positioned[complete.cases(Phen_GRAMMARed_Ctrl_positioned),]

# For trans-eQTL 
# We should not exclude the mitochondrial genes - there might be distal genetic variation
# These (8728:8763) have not been annotated in the past and are NA
# The mapping software will spit out errors if this is NA
# Since the variants are all trans, the positions of the mitochondrial genes don't matter and are set to dummy values
Phen_GRAMMARed_Ctrl_trans_eQTL <- Phen_GRAMMARed_Ctrl_positioned
View(Phen_GRAMMARed_Ctrl_trans_eQTL) 
Phen_GRAMMARed_Ctrl_trans_eQTL[8728:8763,1]  <- "5Mit"
Phen_GRAMMARed_Ctrl_trans_eQTL[8728:8763,"start"] <- 8728:8763 
Phen_GRAMMARed_Ctrl_trans_eQTL[8728:8763,"end"] <- 8729:8764 

# For trans-veQTL
# We need to make dummy positions and chromosomes and manually remove the cis-SNPs after getting back the results
# This has been done for the SNP files as well
Phen_GRAMMARed_Ctrl_trans_veQTL <- Phen_GRAMMARed_Ctrl_positioned
Phen_GRAMMARed_Ctrl_trans_veQTL$`#chr` <- '2L'
Phen_GRAMMARed_Ctrl_trans_veQTL$start <- 1
Phen_GRAMMARed_Ctrl_trans_veQTL$end <- 1
View(Phen_GRAMMARed_Ctrl_trans_veQTL)

#####################################
##### (2b) HS population ############
#####################################

# Same rationale as before
Phen_GRAMMARed_HS_positioned <- merge(Gene_body_locations,Phen_GRAMMARed_HS,all.y = T,by = 'Gene')
Phen_GRAMMARed_HS_positioned[,1:4] <- Phen_GRAMMARed_HS_positioned[,c(4,2,3,1)]
colnames(Phen_GRAMMARed_HS_positioned)[1:4] <- c('chr','start','end','phenotype_id')
Phen_GRAMMARed_HS_positioned <- Phen_GRAMMARed_HS_positioned[order(Phen_GRAMMARed_HS_positioned$chr,Phen_GRAMMARed_HS_positioned$start),]
colnames(Phen_GRAMMARed_HS_positioned)[1] <- '#chr'

# For cis-eQTL and -veQTL
Phen_GRAMMARed_HS_cis <- Phen_GRAMMARed_HS_positioned[complete.cases(Phen_GRAMMARed_HS_positioned),]

# For trans-eQTL
Phen_GRAMMARed_HS_trans_eQTL <- Phen_GRAMMARed_HS_positioned
View(Phen_GRAMMARed_HS_trans_eQTL) # 8728:8763 have NA
Phen_GRAMMARed_HS_trans_eQTL[8728:8763,1] <- "5Mit"
Phen_GRAMMARed_HS_trans_eQTL[8728:8763,"start"] <- 8728:8763 
Phen_GRAMMARed_HS_trans_eQTL[8728:8763,"end"] <- 8729:8764 

# For trans-veQTL 
Phen_GRAMMARed_HS_trans_veQTL <- Phen_GRAMMARed_HS_positioned
Phen_GRAMMARed_HS_trans_veQTL$`#chr` <- '2L'
Phen_GRAMMARed_HS_trans_veQTL$start <- 1
Phen_GRAMMARed_HS_trans_veQTL$end <- 1

#############################################################
##### (3) Save the datasets into .bed files #################
#############################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")

write.table(Phen_GRAMMARed_Ctrl_cis,"Ctrl_phenos_cis_mapping.bed",sep='\t',quote=F,row.names=F)

write.table(Phen_GRAMMARed_HS_cis,"HS_phenos_cis_mapping.bed",sep='\t',quote=F,row.names=F)

write.table(Phen_GRAMMARed_Ctrl_trans_eQTL,"Ctrl_phenos_trans_eQTL_mapping.bed",sep='\t',quote=F,row.names=F)

write.table(Phen_GRAMMARed_HS_trans_eQTL,"HS_phenos_trans_eQTL_mapping.bed",sep='\t',quote=F,row.names=F)

write.table(Phen_GRAMMARed_Ctrl_trans_veQTL,"Ctrl_phenos_trans_veQTL_mapping.bed",sep='\t',quote=F,row.names=F)

write.table(Phen_GRAMMARed_HS_trans_veQTL,"HS_phenos_trans_veQTL_mapping.bed",sep='\t',quote=F,row.names=F)
