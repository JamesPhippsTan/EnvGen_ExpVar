# Preparing phenotype files for use in the mapping with GRAMMAR + veQTL mapper (GraVe)

# Applying the voom transformation and regressing out covariates

# Last updated: 7/11/24

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(limma)
library(tidyr)
library(dplyr)
library(edgeR)
library(tibble)

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
load(file='1_Phenotypes_Covariate_Corrected_GraVe.RData')

#####################
##### Datasets ######
#####################

# Phenotypes - Expression 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
raw.counts <- read.table("RawCounts_noY_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt",h=T,check.names = F)

# Covariates - Environment and batch effects 
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)
cov <- tibble::column_to_rownames(head.info, "id")
cov2 <- cov[c(colnames(raw.counts)),] 

# The TMM and library size normalized data carried out by voom() will be used for the log(variance) mapping by hlmm()
# Thus, the following significant surrogate variables will also be included in the covariate matrix
tmm.voom.sv1 <- read.table('TMM_Voom_sv1_10.txt')
tmm.voom.sv1 <- tmm.voom.sv1$x
tmm.voom.sv2 <- read.table('TMM_Voom_sv2_10.txt')
tmm.voom.sv2 <- tmm.voom.sv2$x
tmm.voom.sv3 <- read.table('TMM_Voom_sv3_10.txt')
tmm.voom.sv3 <- tmm.voom.sv3$x
tmm.voom.sv4 <- read.table('TMM_Voom_sv4_10.txt')
tmm.voom.sv4 <- tmm.voom.sv4$x

# Bimodal genes to be removed
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
bimodal_genes <- read.csv('Bimodal_Genes.csv',header = T)

#################################################################################
##### Phenotype matrix - TMM and library-size normalized counts using voom ######
#################################################################################

# Calculate TMM factors 
raw.counts.list <- DGEList(counts=raw.counts)
TMM.counts <- calcNormFactors(raw.counts.list, method = 'TMM') 

# Preparing independent variable and batch variables
cov2$treatment <- as.factor(cov2$treatment)
cov2$RNAlibBatch <- as.factor(cov2$RNAlibBatch)
cov2$RNAseqBatch <- as.factor(cov2$RNAseqBatch)
cov2$egglayBatch<- as.factor(cov2$egglayBatch)
cov2$platingBatch <- as.factor(cov2$platingBatch)
cov2$well <- as.factor(cov2$well)

treatment <- cov2$treatment
rnalib <- cov2$RNAlibBatch
rnaseq <- cov2$RNAseqBatch
egg <-cov2$egglayBatch
plate <-cov2$platingBatch
well <-cov2$well

contrasts(rnalib) <- contr.sum(levels(rnalib)) # sums to zero contrast
contrasts(rnaseq) <- contr.sum(levels(rnaseq))
contrasts(egg) <- contr.sum(levels(egg))
contrasts(plate) <- contr.sum(levels(plate))
contrasts(well) <- contr.sum(levels(well))

# Design and covariate matrices
TMM.Voom.design <- model.matrix(~rnalib+rnaseq+egg+plate+well+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4+treatment)
TMM.Voom.covariates <- model.matrix(~rnalib+rnaseq+egg+plate+well+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4)

# Transform count data while keeping design
# First log2(cpm+0.5). Then it does voom keeping covariates separated
voom.object <- voom(TMM.counts, design = TMM.Voom.design, plot=T)

# Remove the batch effects
nullmodel <- model.matrix(~treatment)
voom.counts.adj <- limma::removeBatchEffect(voom.object, 
                                            covariates = TMM.Voom.covariates[,-1],  
                                            design = nullmodel) 
voom.counts.adj <- as.data.frame(voom.counts.adj)

# Preparing the final expression matrix 
Expression_Matrix <- as.data.frame(t(voom.counts.adj))
Expression_Matrix_2 <- rownames_to_column(Expression_Matrix)
Expression_Matrix_2$rowname =gsub("_", "",rownames(cov2)) 
# removal of _ is necessary for some mapping steps

# GRAMMAR requires the first column to be id and second column to be sex
Expression_Matrix_3 <- cbind(Expression_Matrix_2$rowname,Expression_Matrix_2)
colnames(Expression_Matrix_3)[1] <- 'id' 
colnames(Expression_Matrix_3)[2] <- 'sex'
Expression_Matrix_3[2] <- 0 # all females so all 0s

# Rearrange
GraVe_phenotypes <- Expression_Matrix_3[order(Expression_Matrix_3[,1]),]
rownames(GraVe_phenotypes) <- 1:nrow(GraVe_phenotypes)

# Remove bimodal genes 
GraVe_phenotypes_bimodalless <- GraVe_phenotypes[,setdiff(colnames(GraVe_phenotypes),bimodal_genes$Gene_name)]

# Get the environments for each individual to subset the vcf and phenotype dataframe 
GraVe_IDs_environments <- data.frame(id=gsub("_", "",rownames(cov2)),environment=cov2$treatment)
# removal of _ is necessary for some mapping steps
ctrl_pop <- subset(GraVe_IDs_environments,environment==1)$id
hs_pop <- subset(GraVe_IDs_environments,environment==2)$id

# Split by control and hs
rownames(GraVe_phenotypes_bimodalless) <- GraVe_phenotypes_bimodalless$id
GraVe_phenotypes_Ctrl <- GraVe_phenotypes_bimodalless[ctrl_pop,]
GraVe_phenotypes_HS <- GraVe_phenotypes_bimodalless[hs_pop,]

# Gene - phenotype number mapping for confirmation purposes
GraVe_phenotypes_map <- data.frame(c(1:length(colnames(Expression_Matrix))),colnames(Expression_Matrix))
colnames(GraVe_phenotypes_map) <- c('Phenotype_Number','Gene_ID')

##################################
##### Save the final matrix ######
##################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")

# Get the phenotypes for GRAMMAR
write.table(GraVe_phenotypes_Ctrl,'GraVe_phenotypes_Ctrl.txt',quote=F,row.names=F)
write.table(GraVe_phenotypes_HS,'GraVe_phenotypes_HS.txt',quote=F,row.names=F)
write.table(GraVe_phenotypes_map,'GraVe_phenotypes_map.txt',quote=F,row.names=F)

# Get the ids and environments for easy subsetting
write.table(GraVe_IDs_environments,'GraVe_IDs_environments.txt',quote=F,row.names=F)
write.table(ctrl_pop,'Ctrl_samples.txt',quote=F,row.names=F,col.names = F)
write.table(hs_pop,'HS_samples.txt',quote=F,row.names=F,col.names = F)

save.image(file='1_Phenotypes_Covariate_Corrected_GraVe.RData')
