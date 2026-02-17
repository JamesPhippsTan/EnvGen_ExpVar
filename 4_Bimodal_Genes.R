# Transcriptome variance analyses of Adult Netherlands D.mel heads on normal vs high sugar diets
# Bimodality checks for different batch effect corrections
# MAD as variability metric

# Last Updated: 7/8/24


#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(DESeq2)
library(dplyr)
library(edgeR)
library(limma)
library(tidyr)
library(car)
library(matrixStats)
library(diptest)
library(tibble)
library(scatterplot3d)

# Working Directory 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")

# Load saved script environment
#load(file='Bimodality_Batches.RData')


#####################
##### Datasets ######
#####################

# Metadata on all samples, e.g., conditions and batches
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)

# Surrogate variables, estimated in the NvsHS_SurrogateVariables R scripts
tmm.voom.sv1 <- read.table('TMM_Voom_sv1_10.txt')
tmm.voom.sv1 <- tmm.voom.sv1$x
tmm.voom.sv2 <- read.table('TMM_Voom_sv2_10.txt')
tmm.voom.sv2 <- tmm.voom.sv2$x
tmm.voom.sv3 <- read.table('TMM_Voom_sv3_10.txt')
tmm.voom.sv3 <- tmm.voom.sv3$x
tmm.voom.sv4 <- read.table('TMM_Voom_sv4_10.txt')
tmm.voom.sv4 <- tmm.voom.sv4$x

vst.sv1 <- read.table('VST_sv1_10.txt')
vst.sv1 <- vst.sv1$x
vst.sv2 <- read.table('VST_sv2_10.txt')
vst.sv2 <- vst.sv2$x
vst.sv3 <- read.table('VST_sv3_10.txt')
vst.sv3 <- vst.sv3$x

# Based on PC inspection - Include 3 SVs for VST; 4 SVs for voom

# Expression - raw counts of filtered samples (see MakingGeneExpressionMatrix_head_HS&CTRL.R for filtering criteria)
raw.counts <- read.table("RawCounts_noY_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt",h=T,check.names = F)

# Only metadata on samples with counts
cov <- tibble::column_to_rownames(head.info, "id")
cov2 <- cov[c(colnames(raw.counts)),]

# Merge expression and metadata 
raw.counts.t <- data.frame(t(raw.counts)) 
raw.counts.t <- tibble::rownames_to_column(raw.counts.t, "id")
head.data <- merge(raw.counts.t,head.info,by = 'id')
head.data <- tibble::column_to_rownames(head.data, "id")


################################
##### (1) Transformations ######
################################

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

# nullmodel - the biological factor of interest
nullmodel <- model.matrix(~treatment)

# Defining indices to subset the total dataset into Ctrl and HS downstream
conditions <- factor(t(head.data$treatment))
conditionsLevel <- levels(conditions)

############################
##### TMM Voom #############
############################

# Calculate TMM factors 
raw.counts.list <- DGEList(counts=raw.counts)
TMM.counts <- calcNormFactors(raw.counts.list, method = 'TMM') 

# Transform count data while keeping design
# First log2(cpm+0.5). Then it does voom keeping covariates separated
TMM.Voom.design <- model.matrix(~rnalib+rnaseq+egg+plate+well+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4+treatment)
voom.object <- voom(TMM.counts, design = TMM.Voom.design, plot=T)

# For tests that do not accept batch effects as covariates:
# Correct counts for weights and batch effects but not for treatment
TMM.Voom.covariates <- model.matrix(~rnalib+rnaseq+egg+plate+well+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4)
voom.counts.bc <- limma::removeBatchEffect(voom.object, 
                                           covariates = TMM.Voom.covariates[,-1],  
                                           design = nullmodel) 
voom.counts.bc <- as.data.frame(voom.counts.bc)

# Split data by control or high sugar
voomdataN <- voom.counts.bc[,c(which(conditions==conditionsLevel[1]))]
voomdataHS <- voom.counts.bc[,c(which(conditions==conditionsLevel[2]))]

###################
##### VST #########
###################

# Apply variance-stabilizing transformation
data <- DESeqDataSetFromMatrix(countData = raw.counts,
                               colData = cov2,
                               design =  ~treatment) 
vst.data = vst(data, blind=F) 

# For tests that do not accept batch effects as covariates:
# Correct counts for batch effects but not for treatment
vst.counts <- assay(vst.data)

# Remove batch effects 
VST.covariates <- model.matrix(~rnalib+rnaseq+egg+plate+well+vst.sv1+vst.sv2+vst.sv3)
vst.counts.bc  <- limma::removeBatchEffect(vst.counts, covariates = VST.covariates[,-1],  design = nullmodel) 
vst.counts.bc.t <- t(vst.counts.bc)

# Split data by control or high sugar
VSTdataN <- vst.counts.bc[,c(which(conditions==conditionsLevel[1]))]
VSTdataHS <- vst.counts.bc[,c(which(conditions==conditionsLevel[2]))]

#############################################################################
##### (2) Checking for bimodality in count distributions within conditions ##
#############################################################################

############################################################################
##### (2a) After all batch effect and surrogate variable correction ########
############################################################################

# It is suspected that there maybe bimodal distributions within conditions due to experimental procedures
# For example, the head separation procedure may have removed antennae from samples and not others
# As such, these may be identified using a Hartigan's Dip test and be excluded at some point during downstream analyses

# Prepare dataframe
Diptest_pvals <- data.frame(Gene_name = NA*8779,
                            Ctrl_pval_voomdata = NA*8779,
                            HS_pval_voomdata = NA*8779,
                            Ctrl_pval_vstdata = NA*8779,
                            HS_pval_vstdata = NA*8779)

# Transpose count data 
voomdataN.t <- t(voomdataN)
voomdataHS.t <- t(voomdataHS)
VSTdataN.t <- t(VSTdataN)
VSTdataHS.t <- t(VSTdataHS)

# Use Hartigan's dip test to ascertain if there are bimodally-distributed genes
for (gene_index in 1:8779){
  gene_name <- rownames(VSTdataN)[gene_index]
  Diptest_pvals[gene_index,'Gene_name'] <- gene_name
  Diptest_pvals[gene_index,'Ctrl_pval_voomdata'] <- dip.test(voomdataN.t[,gene_name])$p.value
  Diptest_pvals[gene_index,'HS_pval_voomdata'] <- dip.test(voomdataHS.t[,gene_name])$p.value
  Diptest_pvals[gene_index,'Ctrl_pval_vstdata'] <- dip.test(VSTdataN.t[,gene_name])$p.value
  Diptest_pvals[gene_index,'HS_pval_vstdata'] <- dip.test(VSTdataHS.t[,gene_name])$p.value
}

# FDR correction
Diptest_pvals[,'Ctrl_pval_voomdata'] <- p.adjust(Diptest_pvals[,'Ctrl_pval_voomdata'],method='BH')
Diptest_pvals[,'HS_pval_voomdata'] <- p.adjust(Diptest_pvals[,'HS_pval_voomdata'],method='BH')
Diptest_pvals[,'Ctrl_pval_vstdata'] <- p.adjust(Diptest_pvals[,'Ctrl_pval_vstdata'],method='BH')
Diptest_pvals[,'HS_pval_vstdata'] <- p.adjust(Diptest_pvals[,'HS_pval_vstdata'],method='BH')

# Getting those that are significantly bimodal (less than 0.05)
Bimodal_genes <- subset(Diptest_pvals,Ctrl_pval_vstdata < 0.05 | HS_pval_vstdata < 0.05 | Ctrl_pval_voomdata < 0.05 | HS_pval_voomdata < 0.05)
nrow(Bimodal_genes) # 16 genes

# Dissecting these by category and transformation
Bimodal_genes_voomdata <- subset(Diptest_pvals,Ctrl_pval_voomdata < 0.05 | HS_pval_voomdata < 0.05)
nrow(Bimodal_genes_voomdata)
# There are 12 genes in the voom data 
voom_ctrl <- subset(Bimodal_genes,Ctrl_pval_voomdata < 0.05) # 12
voom_hs <- subset(Bimodal_genes,HS_pval_voomdata < 0.05) # 2
intersect(voom_ctrl$Gene_name,voom_hs$Gene_name)
# All genes intersect

Bimodal_genes_vstdata <- subset(Diptest_pvals,Ctrl_pval_vstdata < 0.05 | HS_pval_vstdata < 0.05)
nrow(Bimodal_genes_vstdata)
# There are 14 genes in the vst data 
vst_ctrl <- subset(Bimodal_genes,Ctrl_pval_vstdata < 0.05) # 14
vst_hs <- subset(Bimodal_genes,HS_pval_vstdata < 0.05) # 1
intersect(vst_ctrl$Gene_name,vst_hs$Gene_name) # all genes intersect

# intersect between conditions 
intersect(vst_ctrl$Gene_name,voom_ctrl$Gene_name) # 10
intersect(vst_hs$Gene_name,voom_hs$Gene_name) # none
length(intersect(Bimodal_genes_voomdata$Gene_name,Bimodal_genes_vstdata$Gene_name))
# 10 of these genes overlap 

# Plot these out to visualize if the test did its job
par(mfrow = c(2, 2))
for (number in 1:nrow(voom_ctrl)) {
  print(nrow(voom_ctrl))
  bimodal_gene <- voom_ctrl$Gene_name[number]
  plot(density(voomdataN.t[,bimodal_gene]),main = paste0('Ctrl_Voom_',bimodal_gene))
  plot(density(voomdataHS.t[,bimodal_gene]),main = paste0('HS_Voom_',bimodal_gene))
  plot(density(VSTdataN.t[,bimodal_gene]),main = paste0('Ctrl_VST_',bimodal_gene))
  plot(density(VSTdataHS.t[,bimodal_gene]),main = paste0('HS_VST_',bimodal_gene))
}

for (number in 1:nrow(voom_hs)) {
  print(nrow(voom_hs))
  bimodal_gene <- voom_hs$Gene_name[number]
  plot(density(voomdataN.t[,bimodal_gene]),main = paste0('Ctrl_Voom_',bimodal_gene))
  plot(density(voomdataHS.t[,bimodal_gene]),main = paste0('HS_Voom_',bimodal_gene))
  plot(density(VSTdataN.t[,bimodal_gene]),main = paste0('Ctrl_VST_',bimodal_gene))
  plot(density(VSTdataHS.t[,bimodal_gene]),main = paste0('HS_VST_',bimodal_gene))}

for (number in 1:nrow(vst_ctrl)) {
  print(nrow(vst_ctrl))
  bimodal_gene <- vst_ctrl$Gene_name[number]
  plot(density(voomdataN.t[,bimodal_gene]),main = paste0('Ctrl_Voom_',bimodal_gene))
  plot(density(voomdataHS.t[,bimodal_gene]),main = paste0('HS_Voom_',bimodal_gene))
  plot(density(VSTdataN.t[,bimodal_gene]),main = paste0('Ctrl_VST_',bimodal_gene))
  plot(density(VSTdataHS.t[,bimodal_gene]),main = paste0('HS_VST_',bimodal_gene))}

for (number in 1:nrow(vst_hs)) {
  bimodal_gene <- vst_hs$Gene_name[number]
  plot(density(voomdataN.t[,bimodal_gene]),main = paste0('Ctrl_Voom_',bimodal_gene))
  plot(density(voomdataHS.t[,bimodal_gene]),main = paste0('HS_Voom_',bimodal_gene))
  plot(density(VSTdataN.t[,bimodal_gene]),main = paste0('Ctrl_VST_',bimodal_gene))
  plot(density(VSTdataHS.t[,bimodal_gene]),main = paste0('HS_VST_',bimodal_gene))}
# Yeah these definitely look bimodal
# Note: 10403 only not bimodal in voomdataHS

# Example gene where voom data is bimodal but vst is not 
par(mfrow = c(2, 2))
plot(density(voomdataN.t[,'FBgn0036489']),main = paste0('Ctrl_Voom_','FBgn0036489'))
plot(density(voomdataHS.t[,'FBgn0036489']),main = paste0('HS_Voom_','FBgn0036489'))
plot(density(VSTdataN.t[,'FBgn0036489']),main = paste0('Ctrl_VST_','FBgn0036489'))
plot(density(VSTdataHS.t[,'FBgn0036489']),main = paste0('HS_VST_','FBgn0036489'))

# Examine genes of interest
highveQTLgenes <- c('FBgn0037975',
               'FBgn0036184',
               'FBgn0038733',
               'FBgn0087040',
               'FBgn0037973')
for (highveqtlgene in highveQTLgenes){
par(mfrow = c(2, 2))
plot(density(voomdataN.t[,highveqtlgene]),main = paste0('Ctrl_Voom_',highveqtlgene))
plot(density(voomdataHS.t[,highveqtlgene]),main = paste0('HS_Voom_',highveqtlgene))
plot(density(VSTdataN.t[,highveqtlgene]),main = paste0('Ctrl_VST_',highveqtlgene))
plot(density(VSTdataHS.t[,highveqtlgene]),main = paste0('HS_VST_',highveqtlgene))
}
# Only 2 of the 5 somewhat bimodal, though they exhibit a large range of values

#############################################################
##### (2d) Batch effect correction only #####################
#############################################################

bc_only_incl_well <- model.matrix(~rnalib+rnaseq+egg+plate+well)
vst.counts.bc3  <- limma::removeBatchEffect(vst.counts, covariates = bc_only_incl_well[,-1],  design = nullmodel) 
VSTdataN3 <- t(vst.counts.bc3[,c(which(conditions==conditionsLevel[1]))])
VSTdataHS3 <- t(vst.counts.bc3[,c(which(conditions==conditionsLevel[2]))])

voom.counts.adj3 <- limma::removeBatchEffect(voom.object, 
                                             covariates = bc_only_incl_well[,-1],  
                                             design = nullmodel) 
voomdataN3 <- t(voom.counts.adj3[,c(which(conditions==conditionsLevel[1]))])
voomdataHS3 <- t(voom.counts.adj3[,c(which(conditions==conditionsLevel[2]))])

# Is the bimodality present before surrogate variable correction?
Diptest_pvals.noSVinclWell <- data.frame(Gene_name = NA*8779,
                                         Ctrl_pval_voomdata = NA*8779,
                                         HS_pval_voomdata = NA*8779,
                                         Ctrl_pval_vstdata = NA*8779,
                                         HS_pval_vstdata = NA*8779)

# Use Hartigan's dip test to ascertain if there are bimodally-distributed genes
for (gene_index in 1:8779){
  gene_name <- colnames(VSTdataN3)[gene_index]
  Diptest_pvals.noSVinclWell[gene_index,'Gene_name'] <- gene_name
  Diptest_pvals.noSVinclWell[gene_index,'Ctrl_pval_voomdata'] <- dip.test(voomdataN3[,gene_name])$p.value
  Diptest_pvals.noSVinclWell[gene_index,'HS_pval_voomdata'] <- dip.test(voomdataHS3[,gene_name])$p.value
  Diptest_pvals.noSVinclWell[gene_index,'Ctrl_pval_vstdata'] <- dip.test(VSTdataN3[,gene_name])$p.value
  Diptest_pvals.noSVinclWell[gene_index,'HS_pval_vstdata'] <- dip.test(VSTdataHS3[,gene_name])$p.value
}

# View
Diptest_pvals.noSVinclWell[,'Ctrl_pval_voomdata'] <- p.adjust(Diptest_pvals.noSVinclWell[,'Ctrl_pval_voomdata'],method='BH')
Diptest_pvals.noSVinclWell[,'HS_pval_voomdata'] <- p.adjust(Diptest_pvals.noSVinclWell[,'HS_pval_voomdata'],method='BH')
Diptest_pvals.noSVinclWell[,'Ctrl_pval_vstdata'] <- p.adjust(Diptest_pvals.noSVinclWell[,'Ctrl_pval_vstdata'],method='BH')
Diptest_pvals.noSVinclWell[,'HS_pval_vstdata'] <- p.adjust(Diptest_pvals.noSVinclWell[,'HS_pval_vstdata'],method='BH')


# Getting those that are significantly bimodal (less than 0.05)
Bimodal_genes.noSVinclWell_sig <- subset(Diptest_pvals.noSVinclWell,Ctrl_pval_vstdata < 0.1 | HS_pval_vstdata < 0.1 | Ctrl_pval_voomdata < 0.1 | HS_pval_voomdata < 0.1)
View(Bimodal_genes.noSVinclWell_sig)
nrow(Bimodal_genes.noSVinclWell_sig)
# 86 genes without FDR

# Dissecting these by category and transformation
Bimodal_genes_voomdata.noSVinclWell <- subset(Diptest_pvals.noSVinclWell,Ctrl_pval_voomdata < 0.1 | HS_pval_voomdata < 0.1)
nrow(Bimodal_genes_voomdata.noSVinclWell)
# 31 genes
Bimodal_genes_vstdata.noSVinclWell <- subset(Diptest_pvals.noSVinclWell,Ctrl_pval_vstdata < 0.1 | HS_pval_vstdata < 0.1)
nrow(Bimodal_genes_vstdata.noSVinclWell)
# 74 genes
length(intersect(Bimodal_genes_voomdata.noSVinclWell$Gene_name,Bimodal_genes_vstdata.noSVinclWell$Gene_name))
# 19 genes 

# Plotting these out
voom_unique_genes <- setdiff(Bimodal_genes_voomdata.noSVinclWell$Gene_name,Bimodal_genes_vstdata.noSVinclWell$Gene_name)
length(voom_unique_genes) #12
for (index in 1:length(voom_unique_genes)) {
  bimodal_gene <- voom_unique_genes[index]
  par(mfrow=c(2,2))
  plot(density(voomdataN3[,bimodal_gene]),main = paste0('Ctrl_Voom_',bimodal_gene))
  plot(density(voomdataHS3[,bimodal_gene]),main = paste0('HS_Voom_',bimodal_gene))
  plot(density(VSTdataN3[,bimodal_gene]),main = paste0('Ctrl_VST_',bimodal_gene))
  plot(density(VSTdataHS3[,bimodal_gene]),main = paste0('HS_VST_',bimodal_gene))
}

vst_unique_genes <- setdiff(Bimodal_genes_vstdata.noSVinclWell$Gene_name,Bimodal_genes_voomdata.noSVinclWell$Gene_name)
length(vst_unique_genes) #55
for (bimodal_gene in vst_unique_genes) {
  par(mfrow=c(2,2))
  plot(density(voomdataN3[,bimodal_gene]),main = paste0('Ctrl_Voom_',bimodal_gene))
  plot(density(voomdataHS3[,bimodal_gene]),main = paste0('HS_Voom_',bimodal_gene))
  plot(density(VSTdataN3[,bimodal_gene]),main = paste0('Ctrl_VST_',bimodal_gene))
  plot(density(VSTdataHS3[,bimodal_gene]),main = paste0('HS_VST_',bimodal_gene))
}

# Testing coexpression - is the bimodality of these genes due to a common cause?

# Let's start with this gene, which had the highest loading on the residual PC after batch effect correction
PC_1_structure_gene <- "FBgn0051813" # other name is nemure
for (bimodal_gene in Bimodal_genes.noSVinclWell_sig$Gene_name){
  par(mfrow=c(2,2))
  plot(y=voomdataN3[,bimodal_gene],x=voomdataN3[,PC_1_structure_gene],main='voom_control_coexpression',xlab=PC_1_structure_gene,ylab=bimodal_gene)
  plot(y=voomdataHS3[,bimodal_gene],x=voomdataHS3[,PC_1_structure_gene],main='voom_hs_coexpression',xlab=PC_1_structure_gene,ylab=bimodal_gene)
  plot(y=VSTdataN3[,bimodal_gene],x=VSTdataN3[,PC_1_structure_gene],main='vst_control_coexpression',xlab=PC_1_structure_gene,ylab=bimodal_gene)
  plot(y=VSTdataHS3[,bimodal_gene],x=VSTdataHS3[,PC_1_structure_gene],main='vst_hs_coexpression',xlab=PC_1_structure_gene,ylab=bimodal_gene)
}
# Some are coexpressed, some have subtle splits, others have two independent sources of bimodality! (4-way split pattern)

# One of the 4-way-split genes is methuselah8. 
# It still remains bimodal after the 4-sv correction and hence why I wasn't convinced it was just antennae being missing.
# Let's see how it correlates with the other genes in the list
methuselah8 <- "FBgn0052475"
for (bimodal_gene in Bimodal_genes.noSVinclWell_sig$Gene_name){
  par(mfrow=c(2,2))
  plot(y=voomdataN3[,bimodal_gene],x=voomdataN3[,methuselah8],main='voom_control_coexpression',xlab=methuselah8,ylab=bimodal_gene)
  plot(y=voomdataHS3[,bimodal_gene],x=voomdataHS3[,methuselah8],main='voom_hs_coexpression',xlab=methuselah8,ylab=bimodal_gene)
  plot(y=VSTdataN3[,bimodal_gene],x=VSTdataN3[,methuselah8],main='vst_control_coexpression',xlab=methuselah8,ylab=bimodal_gene)
  plot(y=VSTdataHS3[,bimodal_gene],x=VSTdataHS3[,methuselah8],main='vst_hs_coexpression',xlab=methuselah8,ylab=bimodal_gene)
}

# Another is Obp19a
Obp19a <- "FBgn0031109"
for (bimodal_gene in Bimodal_genes.noSVinclWell_sig$Gene_name){
  par(mfrow=c(2,2))
  plot(y=voomdataN3[,bimodal_gene],x=voomdataN3[,Obp19a],main='voom_control_coexpression',xlab=Obp19a,ylab=bimodal_gene)
  plot(y=voomdataHS3[,bimodal_gene],x=voomdataHS3[,Obp19a],main='voom_hs_coexpression',xlab=Obp19a,ylab=bimodal_gene)
  plot(y=VSTdataN3[,bimodal_gene],x=VSTdataN3[,Obp19a],main='vst_control_coexpression',xlab=Obp19a,ylab=bimodal_gene)
  plot(y=VSTdataHS3[,bimodal_gene],x=VSTdataHS3[,Obp19a],main='vst_hs_coexpression',xlab=Obp19a,ylab=bimodal_gene)
}

# Ribosomal protein 49 as a negative control
plot(y=voomdataN3[,"FBgn0002626"],x=voomdataN3[,Obp19a],main='voom_control_coexpression',xlab=Obp19a,ylab=bimodal_gene)
plot(y=voomdataN3[,"FBgn0002626"],x=voomdataN3[,PC_1_structure_gene],main='voom_control_coexpression',xlab=Obp19a,ylab=bimodal_gene)
plot(y=voomdataN3[,PC_1_structure_gene],x=voomdataN3[,Obp19a],main='voom_control_coexpression',xlab=Obp19a,ylab=bimodal_gene)

# How many clusters?
par(mfrow=c(1,1))
pca.bimodals <- prcomp(t(voom.counts.adj3[Bimodal_genes.noSVinclWell_sig$Gene_name,]))
plot(pca.bimodals$x[,1],pca.bimodals$x[,2])

pca.w = NULL
for (i in 1:length(pca.bimodals$sdev)) {
  pca.w[i]= (pca.bimodals$sdev[i])^2/sum(pca.bimodals$sdev^2) 
}
plot(pca.w[1:50])

plot(pca.bimodals$x,pca.bimodals$x[,2])

pc_loadings <- pca.bimodals$rotation
View(pc_loadings)
hist(pc_loadings[,'PC1'])
hist(pc_loadings[,'PC2'])

# Checking that they do not correlate
plot(voomdataN3[,"FBgn0033721"],voomdataN3[,"FBgn0051813"])
plot(voomdataN3[,"FBgn0033721"],voomdataN3[,"FBgn0037414"])
plot(voomdataN3[,"FBgn0033721"],voomdataN3[,"FBgn0037424"])
plot(voomdataN3[,"FBgn0033721"],voomdataN3[,"FBgn0004228"])

###########################################
##### (2d) No batch effect correction #####
###########################################

# Is the bimodality present before batch effect correction?
Diptest_pvals.nocor <- data.frame(Gene_name = NA*8779,
                                  Ctrl_pval_voomdata = NA*8779,
                                  HS_pval_voomdata = NA*8779,
                                  Ctrl_pval_vstdata = NA*8779,
                                  HS_pval_vstdata = NA*8779)
# Transpose count data 
voomdataN.nocor <- t(voom.object$E[,c(which(conditions==conditionsLevel[1]))])
voomdataHS.nocor <- t(voom.object$E[,c(which(conditions==conditionsLevel[2]))])
VSTdataN.nocor <- t(vst.counts[,c(which(conditions==conditionsLevel[1]))])
VSTdataHS.nocor <- t(vst.counts[,c(which(conditions==conditionsLevel[2]))])

# Use Hartigan's dip test to ascertain if there are bimodally-distributed genes
for (gene_index in 1:8779){
  gene_name <- colnames(voomdataN.nocor)[gene_index]
  Diptest_pvals.nocor[gene_index,'Gene_name'] <- gene_name
  Diptest_pvals.nocor[gene_index,'Ctrl_pval_voomdata'] <- dip.test(voomdataN.nocor[,gene_name])$p.value
  Diptest_pvals.nocor[gene_index,'HS_pval_voomdata'] <- dip.test(voomdataHS.nocor[,gene_name])$p.value
  Diptest_pvals.nocor[gene_index,'Ctrl_pval_vstdata'] <- dip.test(VSTdataN.nocor[,gene_name])$p.value
  Diptest_pvals.nocor[gene_index,'HS_pval_vstdata'] <- dip.test(VSTdataHS.nocor[,gene_name])$p.value
}

# FDR correction
Diptest_pvals.nocor[,'Ctrl_pval_voomdata'] <- p.adjust(Diptest_pvals.nocor[,'Ctrl_pval_voomdata'],method='BH')
Diptest_pvals.nocor[,'HS_pval_voomdata'] <- p.adjust(Diptest_pvals.nocor[,'HS_pval_voomdata'],method='BH')
Diptest_pvals.nocor[,'Ctrl_pval_vstdata'] <- p.adjust(Diptest_pvals.nocor[,'Ctrl_pval_vstdata'],method='BH')
Diptest_pvals.nocor[,'HS_pval_vstdata'] <- p.adjust(Diptest_pvals.nocor[,'HS_pval_vstdata'],method='BH')

# Getting those that are significantly bimodal (less than 0.05) before correction
Bimodal_genes.nocor <- subset(Diptest_pvals.nocor,Ctrl_pval_vstdata < 0.05 | HS_pval_vstdata < 0.05 | Ctrl_pval_voomdata < 0.05 | HS_pval_voomdata < 0.05)
nrow(Bimodal_genes.nocor)
# 3451 columns!!!!!

# Dissecting these by category and transformation
Bimodal_genes_voomdata.nocor <- subset(Diptest_pvals.nocor,Ctrl_pval_voomdata < 0.05 | HS_pval_voomdata < 0.05)
nrow(Bimodal_genes_voomdata.nocor)
# There are 619 genes!!!!
Bimodal_genes_vstdata.nocor <- subset(Diptest_pvals.nocor,Ctrl_pval_vstdata < 0.05 | HS_pval_vstdata < 0.05)
nrow(Bimodal_genes_vstdata.nocor)
length(intersect(Bimodal_genes_voomdata.nocor$Gene_name,Bimodal_genes_vstdata.nocor$Gene_name))

######################################################
##### End) Save significantly bimodal genes ##########
######################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
write.csv(Bimodal_genes, 'Bimodal_Genes.csv',row.names = F) 
Bimodal_genes_0.1 <- subset(Diptest_pvals.noSVinclWell,Ctrl_pval_vstdata < 0.1 | HS_pval_vstdata < 0.1 | Ctrl_pval_voomdata < 0.1 | HS_pval_voomdata < 0.1)
write.csv(Bimodal_genes_0.1$Gene_name, 'Bimodal_Genes_0.1.csv',row.names = F) 
Bimodal_genes_0.05 <- subset(Diptest_pvals.noSVinclWell,Ctrl_pval_vstdata < 0.05 | HS_pval_vstdata < 0.05 | Ctrl_pval_voomdata < 0.05 | HS_pval_voomdata < 0.05)
write.csv(Bimodal_genes_0.05$Gene_name, 'Bimodal_Genes_0.05.csv',row.names = F) 

# Save workspace
save.image(file='Bimodality_Batches.RData') 
