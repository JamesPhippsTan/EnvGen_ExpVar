# Removing Y-chromosome genes
# Transcriptome datasets of Adult Netherlands D.mel heads on normal vs high sugar diets

# Script by James Tan
# Adapted from Script 'MakingGeneExpressionMatrix_head_HS&CTRL.R' and 'ExploringBatchEffects' by Dr. Luisa Pallares 
# Last Updated: 6/6/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

# Working Directory 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")


#######################
##### Datasets ########
#######################

# Metadata on all samples, e.g., conditions and batches
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)

# Expression - raw counts of filtered samples (see MakingGeneExpressionMatrix_head_HS&CTRL.R for filtering criteria)
raw.counts <- read.table("RawCounts_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt",h=T,check.names = F)
raw.countsnoY <- read.table("RawCounts_noY_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt",h=T,check.names = F)

################################################################
##### Y chromosome genes - checking expression and removal #####
################################################################

# Only metadata on samples with counts
cov <- tibble::column_to_rownames(head.info, "id")
cov2 <- cov[c(colnames(raw.counts)),]

# Merge expression and metadata 
raw.counts.t <- data.frame(t(raw.counts)) 
raw.counts.t <- tibble::rownames_to_column(raw.counts.t, "id")
head.data <- merge(raw.counts.t,head.info,by = 'id')
head.data <- tibble::column_to_rownames(head.data, "id")

# Condition subsetting
# Defining indices to subset the total dataset into Ctrl and HS downstream
conditions <- factor(t(head.data$treatment))
conditionsLevel <- levels(conditions)
Ctrl_flies <- c(which(conditions==conditionsLevel[1]))
HS_flies <- c(which(conditions==conditionsLevel[2]))

# Condition
cov2$treatment <- as.factor(cov2$treatment)

# Preparing known batch variables
cov2$well <- as.factor(cov2$well)
cov2$plate <- as.factor(cov2$plate)
cov2$RNAlibBatch <- as.factor(cov2$RNAlibBatch)
cov2$RNAseqBatch <- as.factor(cov2$RNAseqBatch)
cov2$egglayBatch<- as.factor(cov2$egglayBatch)
cov2$platingBatch <- as.factor(cov2$platingBatch)
cov2$eclosionBatch <- as.factor(cov2$eclosionBatch)

treatment <- cov2$treatment
well <- cov2$well
plate_number <- cov2$plate
rnalib <- cov2$RNAlibBatch
rnaseq <- cov2$RNAseqBatch
egg <-cov2$egglayBatch
plate <- cov2$platingBatch
eclosion <- cov2$eclosionBatch

contrasts(well) <- contr.sum(levels(well))
contrasts(plate_number) <- contr.sum(levels(plate_number))
contrasts(rnalib) <- contr.sum(levels(rnalib))
contrasts(rnaseq) <- contr.sum(levels(rnaseq))
contrasts(egg) <- contr.sum(levels(egg))
contrasts(plate) <- contr.sum(levels(plate))
contrasts(eclosion) <- contr.sum(levels(eclosion))

# Check expression of these genes. 
raw.counts.Y <- DGEList(counts=raw.counts)
TMM.counts <- calcNormFactors(raw.counts.Y, method = 'TMM') 
voom.object <- voom(TMM.counts, design = NULL, plot=T)
voom.expression <- voom.object$E
covariates <- model.matrix(~rnalib)
logcpm0.5 <- limma::removeBatchEffect(voom.expression, covariates = covariates[,-1],  design = model.matrix(~treatment)) 

# Identified Y-chromosome genes by ShinyGO's genome mapping tool
Ygenes <- c("FBgn0001313","FBgn0046323","FBgn0046697","FBgn0267432","FBgn0267433","FBgn0267449","FBgn0267592")
for (gene in Ygenes){
  hist(logcpm0.5[gene,],main="Expression of Gene")
  plot(logcpm0.5[gene,],logcpm0.5[gene,],main="Expression of Gene")
}
# All are very low - mismapping is the most likely cause

# Remove the genes
NoYgenes <- setdiff(rownames(raw.counts),Ygenes)
raw.countsnoY <- raw.counts[NoYgenes,]
write.table(raw.countsnoY, "RawCounts_noY_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt")
