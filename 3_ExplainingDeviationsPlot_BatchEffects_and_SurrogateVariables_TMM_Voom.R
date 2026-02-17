# Identifying surrogate variables in dataset across different transformations
# Transcriptome datasets of Adult Netherlands D.mel heads on normal vs high sugar diets
# The goal is to infer hidden technical confounding variables present in each transformed expression dataset
# And then regress them out

# TMM transformation 

# Script by James Tan
# Adapted from Script 'MakingGeneExpressionMatrix_head_HS&CTRL.R' and 'ExploringBatchEffects' by Dr. Luisa Pallares 
# Last Updated: 22/1/25


#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(dplyr)
library(limma)
library(edgeR)
library(DESeq2)
library(sva)
library(ggplot2)
library(car)
library(stats)
library(lsr)
library(tibble)
library(diptest)
library(gprofiler2)
library(stringr)
library(rrvgo)
library(ggvenn)

# Working Directory 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")

# Load saved script environment
load(file='ExplainingDevationsPlot.RData')

# Load gProfiler functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('gProfiler_Functions.R')
source('Variability_Functions.R')


#######################
##### Datasets ########
#######################

# Metadata on all samples, e.g., conditions and batches
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)

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


######################################################
##### (1) Which batch effects are redundant? #########
######################################################

cramerscorrelation<-matrix(ncol=8,nrow=8)
chisqpvalue<-matrix(ncol=8,nrow=8)
c<- cov2[,c(1,2,7:9,12:14)]
for (i in 1:8){
  first<- c[,i]
  for (j in 1:8){
    second<- c[,j]
    chisqpvalue[i,j]<-chisq.test(table(first,second))$p.value
    cramerscorrelation[i,j] <- cramersV(table(first,second))
  }
}
colnames(cramerscorrelation)<- names(c)
rownames(cramerscorrelation) <- names(c)
cc <- as.data.frame(cramerscorrelation)
View(cc)
colnames(chisqpvalue)<- names(c)
rownames(chisqpvalue) <- names(c)
cpval <- as.data.frame(chisqpvalue)
View(cpval)

# Conclusion
# Plate has very high correlation with other more specific factors - remove
# Eclosion and plating batch have >0.8 correlation - remove eclosion batch
# Run anova later to see effect size and significance


###############################
##### (2) Transformation ######
###############################

# Calculate TMM factors 
raw.counts.list <- DGEList(counts=raw.counts)
TMM.counts <- calcNormFactors(raw.counts.list, method = 'TMM') 
library_sizes <- TMM.counts$samples$lib.size*TMM.counts$samples$norm.factors

# Normalize count data
voom.object <- voom(TMM.counts, design = NULL, plot=T)
voom.counts <- voom.object$E
voom.countsm <- as.matrix(voom.counts)

# Null model (only treatment) to separate by condition during the batch effect correction steps
nullmodel<-model.matrix(~treatment)

########################################################
##### (3) Remaining batch effect correction ############
########################################################

# Compare counts before and after correcting for batch effects
# Run PCA and plot PCs to visually confirm that each batch effect has been corrected for 
# If you cannot see unique batch-specific clusters in the PCA plots, then the correction has done its job

#############################
##### (3a) voom counts  #####
#############################

pca.raw <- prcomp(t(voom.countsm))

# Check the weights, or percentage of variance explained per PC
pca.raw.w = NULL
for (i in 1:length(pca.raw$sdev)) {
  pca.raw.w[i]= (pca.raw$sdev[i])^2/sum(pca.raw$sdev^2) 
}
plot(pca.raw.w[1:50]) # First 50

# Plot the PCs for variables.
fiftypcs.raw <- pca.raw$x
par(mfrow=c(1,1))
i = 2
j = 1

#plot treatment 
plot(fiftypcs.raw[,i]~ fiftypcs.raw[,j], col="white", pch=16, main="Treatment")
points(fiftypcs.raw[,i][cov2$treatment=="1"]~ fiftypcs.raw[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs.raw[,i][cov2$treatment=="2"]~ fiftypcs.raw[,j][cov2$treatment=='2'],col="red", pch=1)
# Huge PC1 and PC2 structure

# Can't really see
unique_wells <- unique(well)
for (well_number in 1:length(unique_wells)){
  plot(fiftypcs.raw[,i]~ fiftypcs.raw[,j], col="white", pch=16, main=paste0("Well",well_number))
  points(fiftypcs.raw[,i][cov2$treatment=="1"]~ fiftypcs.raw[,j][cov2$treatment=='1'], col="black", pch=1)
  points(fiftypcs.raw[,i][cov2$treatment=="2"]~ fiftypcs.raw[,j][cov2$treatment=='2'],col=rgb(64, 64, 64, maxColorValue = 255), pch=1)
  points(fiftypcs.raw[,i][cov2$well==unique_wells[well_number]]~ fiftypcs.raw[,j][cov2$well==unique_wells[well_number]], col="red", pch=1)
}

#plot sequencing batch
plot(fiftypcs.raw[,i] ~ fiftypcs.raw[,j], col="white", pch=16,  main="Seq batch")
points(fiftypcs.raw[,i][cov2$RNAseqBatch==1]~ fiftypcs.raw[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(fiftypcs.raw[,i][cov2$RNAseqBatch==2]~ fiftypcs.raw[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(fiftypcs.raw[,i][cov2$RNAseqBatch==3]~ fiftypcs.raw[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(fiftypcs.raw[,i][cov2$RNAseqBatch==4]~ fiftypcs.raw[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(fiftypcs.raw[,i][cov2$RNAseqBatch==5]~ fiftypcs.raw[,j][cov2$RNAseqBatch==5],col="grey", pch=1)
# PC1 structure present in batches 2,4,5 but not in batch 1 and 3
# ALL PC2 structure present

#plot library prep batch 
plot(fiftypcs.raw[,i] ~ fiftypcs.raw[,j], col="white", pch=16,  main="Lib prep batch")
points(fiftypcs.raw[,i][cov2$RNAlibBatch==1]~ fiftypcs.raw[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(fiftypcs.raw[,i][cov2$RNAlibBatch==2]~ fiftypcs.raw[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(fiftypcs.raw[,i][cov2$RNAlibBatch==3]~ fiftypcs.raw[,j][cov2$RNAlibBatch==3],col="black", pch=1)
# PC1 structure present in batches 2,3 but not in batch 1
# ALL PC2 structure present

#plot egglay batch 
plot(fiftypcs.raw[,i] ~ fiftypcs.raw[,j], col="white", pch=16,  main="Egglay batch")
c <- sample(colours(),3)
for (e in 1){
  points(fiftypcs.raw[,i][cov2$egglayBatch==e]~ fiftypcs.raw[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 2){
  points(fiftypcs.raw[,i][cov2$egglayBatch==e]~ fiftypcs.raw[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 3){
  points(fiftypcs.raw[,i][cov2$egglayBatch==e]~ fiftypcs.raw[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
# Very subtle - not much

#plot plating day batch
plot(fiftypcs.raw[,i] ~ fiftypcs.raw[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[1]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[5]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
# 1 batch
points(fiftypcs.raw[,i][cov2$platingBatch==ee[2]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[3]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[4]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[6]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
# 1 batch
points(fiftypcs.raw[,i][cov2$platingBatch==ee[7]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)
# Perfectly explains the PC1 structure -> need to regress out for sure

# pca rotation original
pca.raw_loadings <- pca.raw$rotation
View(pca.raw_loadings)

#####################################################
##### (3b) Counts - batch effects without well ######
#####################################################

covariates_nowell <- model.matrix(~rnalib+rnaseq+egg+plate)
counts.adj_nowell <- limma::removeBatchEffect(voom.countsm, covariates = covariates_nowell[,-1],  design = nullmodel) 
pca_nowell <- prcomp(t(counts.adj_nowell))

View(pca_nowell$rotation)

# Look at the percentage of variances explained
pca_nowell.w = NULL
for (i in 1:length(pca_nowell$sdev)) {
  pca_nowell.w[i]= (pca_nowell$sdev[i])^2/sum(pca_nowell$sdev^2) 
}
plot(pca_nowell.w[1:50])

# Let's look at the 1st and 2nd PCs
fiftypcs_nowell <- pca_nowell$x
par(mfrow=c(1,1))
i = 2
j = 1

#plot treatment 
plot(fiftypcs_nowell[,i]~ fiftypcs_nowell[,j], col="white", pch=16, main="Treatment")
points(fiftypcs_nowell[,i][cov2$treatment=="1"]~ fiftypcs_nowell[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs_nowell[,i][cov2$treatment=="2"]~ fiftypcs_nowell[,j][cov2$treatment=='2'],col="red", pch=1)
# Large PC2 and some big PC1 structure

#plot sequencing batch - is even now, but the data is still structured at the first pc
plot(fiftypcs_nowell[,i] ~ fiftypcs_nowell[,j], col="white", pch=16,  main="Seq batch")
points(fiftypcs_nowell[,i][cov2$RNAseqBatch==1]~ fiftypcs_nowell[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(fiftypcs_nowell[,i][cov2$RNAseqBatch==2]~ fiftypcs_nowell[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(fiftypcs_nowell[,i][cov2$RNAseqBatch==3]~ fiftypcs_nowell[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(fiftypcs_nowell[,i][cov2$RNAseqBatch==4]~ fiftypcs_nowell[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(fiftypcs_nowell[,i][cov2$RNAseqBatch==5]~ fiftypcs_nowell[,j][cov2$RNAseqBatch==5],col="grey", pch=1)

#plot library prep batch - is even now, but the data is still structured at the first pc
plot(fiftypcs_nowell[,i] ~ fiftypcs_nowell[,j], col="white", pch=16,  main="Lib prep batch")
points(fiftypcs_nowell[,i][cov2$RNAlibBatch==1]~ fiftypcs_nowell[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(fiftypcs_nowell[,i][cov2$RNAlibBatch==2]~ fiftypcs_nowell[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(fiftypcs_nowell[,i][cov2$RNAlibBatch==3]~ fiftypcs_nowell[,j][cov2$RNAlibBatch==3],col="black", pch=1)

#plot egglay batch - is even now, but the data is still structured at the first pc
plot(fiftypcs_nowell[,i] ~ fiftypcs_nowell[,j], col="white", pch=16,  main="Egglay batch")
c <- sample(colours(),3)
for (e in 1){
  points(fiftypcs_nowell[,i][cov2$egglayBatch==e]~ fiftypcs_nowell[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 2){
  points(fiftypcs_nowell[,i][cov2$egglayBatch==e]~ fiftypcs_nowell[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 3){
  points(fiftypcs_nowell[,i][cov2$egglayBatch==e]~ fiftypcs_nowell[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}

#plot plating day batch - is even now, but the data is still structured at the first pc
plot(fiftypcs_nowell[,i] ~ fiftypcs_nowell[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(fiftypcs_nowell[,i][cov2$platingBatch==ee[1]]~ fiftypcs_nowell[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(fiftypcs_nowell[,i][cov2$platingBatch==ee[2]]~ fiftypcs_nowell[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(fiftypcs_nowell[,i][cov2$platingBatch==ee[3]]~ fiftypcs_nowell[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(fiftypcs_nowell[,i][cov2$platingBatch==ee[4]]~ fiftypcs_nowell[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(fiftypcs_nowell[,i][cov2$platingBatch==ee[5]]~ fiftypcs_nowell[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
points(fiftypcs_nowell[,i][cov2$platingBatch==ee[6]]~ fiftypcs_nowell[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
points(fiftypcs_nowell[,i][cov2$platingBatch==ee[7]]~ fiftypcs_nowell[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)

# Well hardly changes the distribution but is significant nonetheless
sum(pca.raw.w[1:200]) # 56% of the total variance explained by these PCs
m.raw1 <-lm(pca.raw$x[,1:200] ~ treatment + well + rnalib + rnaseq + egg + plate)
m.raw2 <-lm(pca.raw$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate)
anova(m.raw1,m.raw2)

##########################################################
##### (3c) Counts - without batch effects + well #########
##########################################################

covariates <- model.matrix(~rnalib+rnaseq+egg+plate+well)
counts.adj <- limma::removeBatchEffect(voom.countsm, covariates = covariates[,-1],  design = nullmodel) 

pca <- prcomp(t(counts.adj))

# Look at percent of variance explained per PC - 6% now
pca.w = NULL
for (i in 1:length(pca$sdev)) {
  pca.w[i]= (pca$sdev[i])^2/sum(pca$sdev^2) 
}
plot(pca.w[1:50])

fiftypcs <- pca$x
par(mfrow=c(1,1))

# Let's look at the 1st and 2nd PCs
i = 2
j = 1

#plot treatment 
plot(fiftypcs[,i]~ fiftypcs[,j], col="white", pch=16, main="Batchless Voom Data ~ Treatment",ylab = "PC2",xlab = "PC1")
points(fiftypcs[,i][cov2$treatment=="1"]~ fiftypcs[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs[,i][cov2$treatment=="2"]~ fiftypcs[,j][cov2$treatment=='2'],col="red", pch=1)
# PC1 structure still present
# PC2 structure mostly removed

#plot sequencing batch - is even now, but the data is still structured at the first pc
plot(fiftypcs[,i] ~ fiftypcs[,j], col="white", pch=16,  main="Seq batch")
points(fiftypcs[,i][cov2$RNAseqBatch==1]~ fiftypcs[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(fiftypcs[,i][cov2$RNAseqBatch==2]~ fiftypcs[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(fiftypcs[,i][cov2$RNAseqBatch==3]~ fiftypcs[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(fiftypcs[,i][cov2$RNAseqBatch==4]~ fiftypcs[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(fiftypcs[,i][cov2$RNAseqBatch==5]~ fiftypcs[,j][cov2$RNAseqBatch==5],col="grey", pch=1)

#plot library prep batch - is even now, but the data is still structured at the first pc
plot(fiftypcs[,i] ~ fiftypcs[,j], col="white", pch=16,  main="Lib prep batch")
points(fiftypcs[,i][cov2$RNAlibBatch==1]~ fiftypcs[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(fiftypcs[,i][cov2$RNAlibBatch==2]~ fiftypcs[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(fiftypcs[,i][cov2$RNAlibBatch==3]~ fiftypcs[,j][cov2$RNAlibBatch==3],col="black", pch=1)

#plot egglay batch - is even now, but the data is still structured at the first pc
plot(fiftypcs[,i] ~ fiftypcs[,j], col="white", pch=16,  main="Egglay batch")
c <- sample(colours(),3)
for (e in 1){
  points(fiftypcs[,i][cov2$egglayBatch==e]~ fiftypcs[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 2){
  points(fiftypcs[,i][cov2$egglayBatch==e]~ fiftypcs[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 3){
  points(fiftypcs[,i][cov2$egglayBatch==e]~ fiftypcs[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}

#plot plating day batch - is even now, but the data is still structured at the first pc
plot(fiftypcs[,i] ~ fiftypcs[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(fiftypcs[,i][cov2$platingBatch==ee[1]]~ fiftypcs[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(fiftypcs[,i][cov2$platingBatch==ee[2]]~ fiftypcs[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(fiftypcs[,i][cov2$platingBatch==ee[3]]~ fiftypcs[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(fiftypcs[,i][cov2$platingBatch==ee[4]]~ fiftypcs[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(fiftypcs[,i][cov2$platingBatch==ee[5]]~ fiftypcs[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
points(fiftypcs[,i][cov2$platingBatch==ee[6]]~ fiftypcs[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
points(fiftypcs[,i][cov2$platingBatch==ee[7]]~ fiftypcs[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)

# 4 non-redundant batch effects are present
# Looks like the original PC1 is gone and the original PC2 is now PC1. 
# The variance explained is consistent with this too.

##################################
##### (4) SV identification ######
##################################

# Finding surrogate variables - tried a range from 1 - 10
mod <-  model.matrix(~ rnalib+rnaseq+egg+plate+well+treatment)
mod0 <- model.matrix(~ rnalib+rnaseq+egg+plate+well)
sva.sva4 <- sva(voom.countsm, mod, mod0,n.sv=10)
sv1.4<- as.numeric(sva.sva4$sv[,1])

# Is SV1 similar across different trials
#plot(sv1.1 ~ sv1.2)
#cor.test(sv1.1, sv1.2)
#max(sv1.1-sv1.2)
#plot(sv1.1 ~ sv1.3)
#max(sv1.1-sv1.3)
#cor.test(sv1.1, sv1.3)
#plot(sv1.1 ~ sv1.4)
#max(sv1.1-sv1.4)
#cor.test(sv1.1, sv1.4)
# Very much so

# Collecting other svs to see if they may reduce hidden structure in the data further
sv2.4 <- as.numeric(sva.sva4$sv[,2])
sv3.4 <- as.numeric(sva.sva4$sv[,3])
sv4.4 <- as.numeric(sva.sva4$sv[,4])
sv5.4 <- as.numeric(sva.sva4$sv[,5])
sv6.4 <- as.numeric(sva.sva4$sv[,6])
sv7.4 <- as.numeric(sva.sva4$sv[,7])
sv8.4 <- as.numeric(sva.sva4$sv[,8])

##############################################
##### (3b) Counts - batch effects - 1 SV #####
##############################################

covariates2 <- model.matrix(~rnalib+rnaseq+egg+plate+well+sv1.4)
counts.adj2 <- limma::removeBatchEffect(voom.countsm, covariates = covariates2[,-1],  design = nullmodel) 
pca2 <- prcomp(t(counts.adj2))
plot(pca2$x[,1] ~ pca2$x[,2])
pca2.w = NULL
for (i in 1:length(pca2$sdev)) {
  pca2.w[i]= (pca2$sdev[i])^2/sum(pca2$sdev^2) 
}
plot(pca2.w[1:50])

# Check loadings of PC1
pca2_loadings <- pca2$rotation
View(pca2_loadings)
hist(pca2_loadings[,'PC1'])

fiftypcs2 <- pca2$x[, 1:50]

# Check for residual batch effects
mcounts.adj2 <-lm(pca2$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv2.4 + sv4.4 + sv5.4)
mcounts.adj2_well <-lm(pca2$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + well + sv1.4 + sv2.4 + sv2.4 + sv4.4 + sv5.4)
anova(mcounts.adj2)
anova(mcounts.adj2_well)
anova(mcounts.adj2,mcounts.adj2_well)

# Let's look at the 1st and 2nd PCs
i = 2
j = 1
par(mfrow=c(1,1))
#plot treatment - 
plot(fiftypcs2[,i]~ fiftypcs2[,j], col="white", pch=16, main="Treatment")
points(fiftypcs2[,i][cov2$treatment=="1"]~ fiftypcs2[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs2[,i][cov2$treatment=="2"]~ fiftypcs2[,j][cov2$treatment=='2'],col="red", pch=1)
# PC1 and PC2 structure still present (diagonal gap)

#plot sequencing batch
plot(fiftypcs2[,i] ~ fiftypcs2[,j], col="white", pch=16,  main="Seq batch")
points(fiftypcs2[,i][cov2$RNAseqBatch==1]~ fiftypcs2[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(fiftypcs2[,i][cov2$RNAseqBatch==2]~ fiftypcs2[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(fiftypcs2[,i][cov2$RNAseqBatch==3]~ fiftypcs2[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(fiftypcs2[,i][cov2$RNAseqBatch==4]~ fiftypcs2[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(fiftypcs2[,i][cov2$RNAseqBatch==5]~ fiftypcs2[,j][cov2$RNAseqBatch==5],col="grey", pch=1)

#plot library prep batch
plot(fiftypcs2[,i] ~ fiftypcs2[,j], col="white", pch=16,  main="Lib prep batch")
points(fiftypcs2[,i][cov2$RNAlibBatch==1]~ fiftypcs2[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(fiftypcs2[,i][cov2$RNAlibBatch==2]~ fiftypcs2[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(fiftypcs2[,i][cov2$RNAlibBatch==3]~ fiftypcs2[,j][cov2$RNAlibBatch==3],col="black", pch=1)

#plot egglay batch
plot(fiftypcs2[,i] ~ fiftypcs2[,j], col="white", pch=16,  main="Egglay batch")
c <- sample(colours(),3)
for (e in 1){
  points(fiftypcs2[,i][cov2$egglayBatch==e]~ fiftypcs2[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 2){
  points(fiftypcs2[,i][cov2$egglayBatch==e]~ fiftypcs2[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 3){
  points(fiftypcs2[,i][cov2$egglayBatch==e]~ fiftypcs2[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}

#plot plating day batch
plot(fiftypcs2[,i] ~ fiftypcs2[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(fiftypcs2[,i][cov2$platingBatch==ee[1]]~ fiftypcs2[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(fiftypcs2[,i][cov2$platingBatch==ee[2]]~ fiftypcs2[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(fiftypcs2[,i][cov2$platingBatch==ee[3]]~ fiftypcs2[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(fiftypcs2[,i][cov2$platingBatch==ee[4]]~ fiftypcs2[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(fiftypcs2[,i][cov2$platingBatch==ee[5]]~ fiftypcs2[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
points(fiftypcs2[,i][cov2$platingBatch==ee[6]]~ fiftypcs2[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
points(fiftypcs2[,i][cov2$platingBatch==ee[7]]~ fiftypcs2[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)

# Large PC1 and PC2 separation - a diagonal line through the data

##############################################
##### (3c) Counts - batch effects - 2 SVs ####
##############################################

# Correction for batch effects and 2 surrogate variables
covariates3 <- model.matrix(~rnalib+rnaseq+egg+plate+well+sv1.4+sv2.4)
counts.adj3 <- limma::removeBatchEffect(voom.countsm, covariates = covariates3[,-1],  design = nullmodel) 

pca3 <- prcomp(t(counts.adj3))
plot(pca3$x[,1] ~ pca3$x[,2])
pca3.w = NULL
for (i in 1:length(pca3$sdev)) {
  pca3.w[i]= (pca3$sdev[i])^2/sum(pca3$sdev^2) 
}
plot(pca3.w[1:50])

fiftypcs3 <- pca3$x[, 1:50]
par(mfrow=c(1,1))

# Check for residual batch effects
mcounts.adj3 <-lm(pca3$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4)
mcounts.adj3_well <-lm(pca3$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4 + well)
anova(mcounts.adj3)
anova(mcounts.adj3_well)
anova(mcounts.adj3,mcounts.adj3_well)

# Let's look at the 1st and 2nd PCs
i = 2
j = 1

#plot treatment
View(pca3$rotation)
plot(fiftypcs3[,2]~ fiftypcs3[,1], col="white", pch=16, main=paste0("Batchless Voom Data -2 SVs","\n","~ Treatment for All Genes"),ylab = "PC2",xlab = "PC1")
points(fiftypcs3[,2][cov2$treatment=="1"]~ fiftypcs3[,1][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs3[,2][cov2$treatment=="2"]~ fiftypcs3[,1][cov2$treatment=='2'],col="red", pch=1)
# Some PC2 structure still remaining

#plot sequencing batch
plot(fiftypcs3[,i] ~ fiftypcs3[,j], col="white", pch=16,  main="Seq batch")
points(fiftypcs3[,i][cov2$RNAseqBatch==1]~ fiftypcs3[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(fiftypcs3[,i][cov2$RNAseqBatch==2]~ fiftypcs3[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(fiftypcs3[,i][cov2$RNAseqBatch==3]~ fiftypcs3[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(fiftypcs3[,i][cov2$RNAseqBatch==4]~ fiftypcs3[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(fiftypcs3[,i][cov2$RNAseqBatch==5]~ fiftypcs3[,j][cov2$RNAseqBatch==5],col="grey", pch=1)


#plot library prep batch
plot(fiftypcs3[,i] ~ fiftypcs3[,j], col="white", pch=16,  main="Lib prep batch")
points(fiftypcs3[,i][cov2$RNAlibBatch==1]~ fiftypcs3[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(fiftypcs3[,i][cov2$RNAlibBatch==2]~ fiftypcs3[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(fiftypcs3[,i][cov2$RNAlibBatch==3]~ fiftypcs3[,j][cov2$RNAlibBatch==3],col="black", pch=1)

#plot egglay batch
plot(fiftypcs3[,i] ~ fiftypcs3[,j], col="white", pch=16,  main="Egglay batch")
c <- sample(colours(),3)
for (e in 1){
  points(fiftypcs3[,i][cov2$egglayBatch==e]~ fiftypcs3[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 2){
  points(fiftypcs3[,i][cov2$egglayBatch==e]~ fiftypcs3[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 3){
  points(fiftypcs3[,i][cov2$egglayBatch==e]~ fiftypcs3[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}

#plot plating day batch
plot(fiftypcs3[,i] ~ fiftypcs3[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(fiftypcs3[,i][cov2$platingBatch==ee[1]]~ fiftypcs3[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(fiftypcs3[,i][cov2$platingBatch==ee[2]]~ fiftypcs3[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(fiftypcs3[,i][cov2$platingBatch==ee[3]]~ fiftypcs3[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(fiftypcs3[,i][cov2$platingBatch==ee[4]]~ fiftypcs3[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(fiftypcs3[,i][cov2$platingBatch==ee[5]]~ fiftypcs3[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
points(fiftypcs3[,i][cov2$platingBatch==ee[6]]~ fiftypcs3[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
points(fiftypcs3[,i][cov2$platingBatch==ee[7]]~ fiftypcs3[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)

# Large PC2 separation still present

##############################################
##### (3d) Counts - batch effects - 3 SVs ####
##############################################

# Correction for batch effects and 3 surrogate variables
covariates4 <- model.matrix(~rnalib+rnaseq+egg+plate+well+sv1.4+sv2.4+sv3.4)
counts.adj4 <- limma::removeBatchEffect(voom.countsm, covariates = covariates4[,-1],  design = nullmodel) 

pca4 <- prcomp(t(counts.adj4))
plot(pca4$x[,1] ~ pca4$x[,2])
pca4.w = NULL
for (i in 1:length(pca4$sdev)) {
  pca4.w[i]= (pca4$sdev[i])^2/sum(pca4$sdev^2) 
}
plot(pca4.w[1:50])
fiftypcs4 <- pca4$x[, 1:50]

# Check for residual batch effects
mcounts.adj4 <-lm(pca4$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4)
mcounts.adj4_well <-lm(pca4$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4 + well)
anova(mcounts.adj4)
anova(mcounts.adj4_well)
anova(mcounts.adj4,mcounts.adj4_well)

par(mfrow=c(1,1))
i=2
j=1

#plot treatment
View(pca4$rotation)
plot(fiftypcs4[,2]~ fiftypcs4[,1], col="white", pch=16, main=paste0("Batchless Voom Data -3 SVs","\n","~ Treatment for All Genes"),ylab = "PC2",xlab = "PC1")
points(fiftypcs4[,i][cov2$treatment=="1"]~ fiftypcs4[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs4[,i][cov2$treatment=="2"]~ fiftypcs4[,j][cov2$treatment=='2'],col="red", pch=1)
# A lot of PC2 structure still remains! 

#plot sequencing batch
plot(fiftypcs4[,i] ~ fiftypcs4[,j], col="white", pch=16,  main="Seq batch")
points(fiftypcs4[,i][cov2$RNAseqBatch==1]~ fiftypcs4[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(fiftypcs4[,i][cov2$RNAseqBatch==2]~ fiftypcs4[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(fiftypcs4[,i][cov2$RNAseqBatch==3]~ fiftypcs4[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(fiftypcs4[,i][cov2$RNAseqBatch==4]~ fiftypcs4[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(fiftypcs4[,i][cov2$RNAseqBatch==5]~ fiftypcs4[,j][cov2$RNAseqBatch==5],col="grey", pch=1)


#plot library prep batch
plot(fiftypcs4[,i] ~ fiftypcs4[,j], col="white", pch=16,  main="Lib prep batch")
points(fiftypcs4[,i][cov2$RNAlibBatch==1]~ fiftypcs4[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(fiftypcs4[,i][cov2$RNAlibBatch==2]~ fiftypcs4[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(fiftypcs4[,i][cov2$RNAlibBatch==3]~ fiftypcs4[,j][cov2$RNAlibBatch==3],col="black", pch=1)

#plot egglay batch
plot(fiftypcs4[,i] ~ fiftypcs4[,j], col="white", pch=16,  main="Egglay batch")
c <- sample(colours(),3)
for (e in 1){
  points(fiftypcs4[,i][cov2$egglayBatch==e]~ fiftypcs4[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 2){
  points(fiftypcs4[,i][cov2$egglayBatch==e]~ fiftypcs4[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 3){
  points(fiftypcs4[,i][cov2$egglayBatch==e]~ fiftypcs4[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}

#plot plating day batch
plot(fiftypcs4[,i] ~ fiftypcs4[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(fiftypcs4[,i][cov2$platingBatch==ee[1]]~ fiftypcs4[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(fiftypcs4[,i][cov2$platingBatch==ee[2]]~ fiftypcs4[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(fiftypcs4[,i][cov2$platingBatch==ee[3]]~ fiftypcs4[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(fiftypcs4[,i][cov2$platingBatch==ee[4]]~ fiftypcs4[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(fiftypcs4[,i][cov2$platingBatch==ee[5]]~ fiftypcs4[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
points(fiftypcs4[,i][cov2$platingBatch==ee[6]]~ fiftypcs4[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
points(fiftypcs4[,i][cov2$platingBatch==ee[7]]~ fiftypcs4[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)

# Large PC2 separation still present

# Check PC loadings on PC1
pc4_loadings <- pca4$rotation
View(pc4_loadings)
hist(pc4_loadings[,'PC1'])

##############################################
##### (3e) Counts - batch effects - 4 SVs ####
##############################################

# Correction for batch effects and 4 surrogate variables
covariates5 <- model.matrix(~rnalib+rnaseq+egg+plate+well+sv1.4+sv2.4+sv3.4+sv4.4)
counts.adj5 <- limma::removeBatchEffect(voom.countsm, covariates = covariates5[,-1],  design = nullmodel) 

pca5 <- prcomp(t(counts.adj5))
plot(pca5$x[,1] ~ pca5$x[,2])
pca5.w = NULL
for (i in 1:length(pca5$sdev)) {
  pca5.w[i]= (pca5$sdev[i])^2/sum(pca5$sdev^2) 
}
plot(pca5.w[1:50])

fiftypcs5 <- pca5$x[, 1:50]
par(mfrow=c(1,1))

# Check for residual batch effects
mcounts.adj5 <-lm(pca5$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4)
mcounts.adj5_well <-lm(pca5$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4 + well)
anova(mcounts.adj5)
anova(mcounts.adj5_well)
anova(mcounts.adj5,mcounts.adj5_well)

# Let's look at the 1st and 2nd PCs
i = 2
j = 1

#plot treatment
plot(fiftypcs5[,2]~ fiftypcs5[,1], col="white", pch=16, main=paste0("Batchless Voom Data -4 SVs","\n","~ Treatment for All Genes"),ylab = "PC2",xlab = "PC1")
points(fiftypcs5[,i][cov2$treatment=="1"]~ fiftypcs5[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs5[,i][cov2$treatment=="2"]~ fiftypcs5[,j][cov2$treatment=='2'],col="red", pch=1)
# PC2 structure largely gone
# Clean separation between treatments

#plot sequencing batch
plot(fiftypcs5[,i] ~ fiftypcs5[,j], col="white", pch=16,  main="Seq batch")
points(fiftypcs5[,i][cov2$RNAseqBatch==1]~ fiftypcs5[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(fiftypcs5[,i][cov2$RNAseqBatch==2]~ fiftypcs5[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(fiftypcs5[,i][cov2$RNAseqBatch==3]~ fiftypcs5[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(fiftypcs5[,i][cov2$RNAseqBatch==4]~ fiftypcs5[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(fiftypcs5[,i][cov2$RNAseqBatch==5]~ fiftypcs5[,j][cov2$RNAseqBatch==5],col="grey", pch=1)


#plot library prep batch
plot(fiftypcs5[,i] ~ fiftypcs5[,j], col="white", pch=16,  main="Lib prep batch")
points(fiftypcs5[,i][cov2$RNAlibBatch==1]~ fiftypcs5[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(fiftypcs5[,i][cov2$RNAlibBatch==2]~ fiftypcs5[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(fiftypcs5[,i][cov2$RNAlibBatch==3]~ fiftypcs5[,j][cov2$RNAlibBatch==3],col="black", pch=1)

#plot egglay batch
plot(fiftypcs5[,i] ~ fiftypcs5[,j], col="white", pch=16,  main="Egglay batch")
c <- sample(colours(),3)
for (e in 1){
  points(fiftypcs5[,i][cov2$egglayBatch==e]~ fiftypcs5[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 2){
  points(fiftypcs5[,i][cov2$egglayBatch==e]~ fiftypcs5[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 3){
  points(fiftypcs5[,i][cov2$egglayBatch==e]~ fiftypcs5[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}

#plot plating day batch
plot(fiftypcs5[,i] ~ fiftypcs5[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(fiftypcs5[,i][cov2$platingBatch==ee[1]]~ fiftypcs5[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(fiftypcs5[,i][cov2$platingBatch==ee[2]]~ fiftypcs5[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(fiftypcs5[,i][cov2$platingBatch==ee[3]]~ fiftypcs5[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(fiftypcs5[,i][cov2$platingBatch==ee[4]]~ fiftypcs5[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(fiftypcs5[,i][cov2$platingBatch==ee[5]]~ fiftypcs5[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
points(fiftypcs5[,i][cov2$platingBatch==ee[6]]~ fiftypcs5[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
points(fiftypcs5[,i][cov2$platingBatch==ee[7]]~ fiftypcs5[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)
# Large PC1 separation removed

# Check PC loadings on PC1
pc5_loadings <- pca5$rotation
View(pc5_loadings)
hist(pc5_loadings[,'PC1'])


###################################################
##### The deviations plot PC - explaining it ######
###################################################

# PCA plot of the absolute deviations from the median matrix per individual
rowADs <- function(mat, na.rm = TRUE) {
  row_medians <- matrixStats::rowMedians(mat, na.rm = na.rm)
  abs(mat - row_medians)
}
# Calculate deviations from condition-specific medians of each gene
ADs.adj5.Ctrl <- rowADs(counts.adj5[,Ctrl_flies])
ADs.adj5.HS <- rowADs(counts.adj5[,HS_flies])
ADs.adj5.merge <- cbind(ADs.adj5.Ctrl,ADs.adj5.HS)
ADs.adj5 <- ADs.adj5.merge[,colnames(counts.adj5)]
ADpca5 <- prcomp(t(ADs.adj5))
plot(ADpca5$x[,1] ~ ADpca5$x[,2])
ADpca5.w = NULL
for (i in 1:length(ADpca5$sdev)) {
  ADpca5.w[i]= (ADpca5$sdev[i])^2/sum(ADpca5$sdev^2) 
}
plot(ADpca5.w[1:50])
ADfiftypcs5 <- ADpca5$x[, 1:50]

# Plot treatment
i=2
j=1
par(mfrow=c(1,2))
plot(ADfiftypcs5[,2]~ ADfiftypcs5[,1], col="white", pch=16, main=paste0("Batchless Voom Data -4 SVs","\n","~ Treatment for All Genes"),ylab = "PC2",xlab = "PC1")
points(ADfiftypcs5[,i][cov2$treatment=="1"]~ ADfiftypcs5[,j][cov2$treatment=='1'], col="blue", pch=1)
points(ADfiftypcs5[,i][cov2$treatment=="2"]~ ADfiftypcs5[,j][cov2$treatment=='2'],col="red", pch=1)
# PC1 structure present in the decomposition of the residual matrix...
plot(fiftypcs5[,2]~ fiftypcs5[,1], col="white", pch=16, main=paste0("Batchless Voom Data -4 SVs","\n","~ Treatment for All Genes"),ylab = "PC2",xlab = "PC1")
points(fiftypcs5[,i][cov2$treatment=="1"]~ fiftypcs5[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs5[,i][cov2$treatment=="2"]~ fiftypcs5[,j][cov2$treatment=='2'],col="red", pch=1)
# not found in the decomposition of the adjusted count matrix...

ADdataframe1 <- data.frame(ID = rownames(ADfiftypcs5), EV1 = ADfiftypcs5[,'PC1'],EV2 = ADfiftypcs5[,'PC2'])
ADdataframe2 <- data.frame(ID = rownames(cov2), Treatment = cov2[,'treatment'])
ADPCplot_df <- merge(ADdataframe1,ADdataframe2,by='ID',all.x=T)
ADPCplot_df$Diet <- ifelse(ADPCplot_df$Treatment==1,'Control','High Sugar')

ADPCplot <- ggplot(data=ADPCplot_df,aes(x=EV1,y=EV2,group=Diet,col=Diet))+geom_point()+ 
  scale_color_manual(values=c('#C6B49F','#DF9F65')) + 
  geom_point() + 
  ylab(paste0('Absolute Deviations PC2 (',100*round(ADpca5.w[2],3),'%)')) + 
  xlab(paste0('Absolute Deviations PC1 (',100*round(ADpca5.w[1],3),'%)')) +
  theme_classic()+  # Clean theme with larger font
  theme(
    legend.position = c(0.2, 0.9),  # Places legend inside plot (x,y between 0-1)
    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) 
ADPCplot
ADPCplot_nolegend <- ADPCplot+ theme(legend.position = "none")
ADPCplot_nolegend
ggsave(filename = 'Transcriptome absolute deviation PCA.png',plot = ADPCplot,height = 3,width = 3)
ggsave(filename = 'Transcriptome absolute deviation PCA nolegend.png',plot = ADPCplot_nolegend,height = 3,width = 3)
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")

# So what's the deal with this structure?

# 1) There is a particular processing batch with many outlier individuals....
#plot sequencing batch
plot(ADfiftypcs5[,i] ~ ADfiftypcs5[,j], col="white", pch=16,  main="Seq batch")
points(ADfiftypcs5[,i][cov2$RNAseqBatch==1]~ ADfiftypcs5[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(ADfiftypcs5[,i][cov2$RNAseqBatch==2]~ ADfiftypcs5[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(ADfiftypcs5[,i][cov2$RNAseqBatch==3]~ ADfiftypcs5[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(ADfiftypcs5[,i][cov2$RNAseqBatch==4]~ ADfiftypcs5[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(ADfiftypcs5[,i][cov2$RNAseqBatch==5]~ ADfiftypcs5[,j][cov2$RNAseqBatch==5],col="grey", pch=1)
#plot library prep batch
plot(ADfiftypcs5[,i] ~ ADfiftypcs5[,j], col="white", pch=16,  main="Lib prep batch")
points(ADfiftypcs5[,i][cov2$RNAlibBatch==1]~ ADfiftypcs5[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(ADfiftypcs5[,i][cov2$RNAlibBatch==2]~ ADfiftypcs5[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(ADfiftypcs5[,i][cov2$RNAlibBatch==3]~ ADfiftypcs5[,j][cov2$RNAlibBatch==3],col="black", pch=1)
#plot egglay batch - there is some batc
plot(ADfiftypcs5[,i] ~ ADfiftypcs5[,j], col="white", pch=16,  main="Egglay batch")
c <- sample(colours(),3)
for (e in 1){
  points(ADfiftypcs5[,i][cov2$egglayBatch==e]~ ADfiftypcs5[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 2){
  points(ADfiftypcs5[,i][cov2$egglayBatch==e]~ ADfiftypcs5[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 3){
  points(ADfiftypcs5[,i][cov2$egglayBatch==e]~ ADfiftypcs5[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
#plot plating day batch
plot(ADfiftypcs5[,i] ~ ADfiftypcs5[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(ADfiftypcs5[,i][cov2$platingBatch==ee[1]]~ ADfiftypcs5[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(ADfiftypcs5[,i][cov2$platingBatch==ee[5]]~ ADfiftypcs5[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
# group 1
points(ADfiftypcs5[,i][cov2$platingBatch==ee[3]]~ ADfiftypcs5[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(ADfiftypcs5[,i][cov2$platingBatch==ee[4]]~ ADfiftypcs5[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(ADfiftypcs5[,i][cov2$platingBatch==ee[2]]~ ADfiftypcs5[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(ADfiftypcs5[,i][cov2$platingBatch==ee[6]]~ ADfiftypcs5[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
# group 2
points(ADfiftypcs5[,i][cov2$platingBatch==ee[7]]~ ADfiftypcs5[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)
# The HS cloud belongs to 1, 5 and 7 exclusively...

# 2) The smallest and largest libraries in this cloud?
# High deviants can be gained from both very high and very low counts
par(mfcol=c(1,1))
plot(ADfiftypcs5[,i] ~ ADfiftypcs5[,j], col="black", pch=16,  main="Plating day batch")
plot(ADfiftypcs5[,1],library_sizes) 
plot(ADfiftypcs5[,1],library_sizes,ylim=c(0,1000000)) 
plot(ADfiftypcs5[,1],library_sizes,ylim=c(1000000,max(library_sizes))) 
# Not really....almost all above 1 million reads though not the higher end of the spectrum
# They are medium-low library size

# 3) Evidence for circadian genes being here? 
# This would provide evidence that such deviants are hig 
ADpca5_loadings <- ADpca5$rotation
View(ADpca5_loadings)
# Some mitochondrial some ribosomal with high positive loadings
# Not different timings on those days

# 4) Created during the batch effect removal procedure?
# What if not including the SVs?
ADs.adj.Ctrl <- rowADs(counts.adj[,Ctrl_flies])
ADs.adj.HS <- rowADs(counts.adj[,HS_flies])
ADs.adj.merge <- cbind(ADs.adj.Ctrl,ADs.adj.HS)
ADs.adj <- ADs.adj.merge[,colnames(counts.adj)]
ADpca <- prcomp(t(ADs.adj))
plot(ADpca$x[,1] ~ ADpca$x[,2])
ADpca.w = NULL
for (i in 1:length(ADpca$sdev)) {
  ADpca.w[i]= (ADpca$sdev[i])^2/sum(ADpca$sdev^2) 
}
plot(ADpca.w[1:0])
ADfiftypcs <- ADpca$x[, 1:50]
par(mfrow=c(1,2))
i = 2
j = 1
plot(ADfiftypcs[,2]~ ADfiftypcs[,1], col="white", pch=16, main='Batch effect correction only',ylab = "Abs Deviation PC2",xlab = "Abs Deviation PC1")
points(ADfiftypcs[,i][cov2$treatment=="1"]~ ADfiftypcs[,j][cov2$treatment=='1'], col="blue", pch=1)
points(ADfiftypcs[,i][cov2$treatment=="2"]~ ADfiftypcs[,j][cov2$treatment=='2'],col="red", pch=1)
# Removing 1 SV still results in the PC1 structure
plot(fiftypcs[,2]~ fiftypcs[,1], col="white", pch=16, main='Batch correction only',ylab = "PC2",xlab = "PC1")
points(fiftypcs[,i][cov2$treatment=="1"]~ fiftypcs[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs[,i][cov2$treatment=="2"]~ fiftypcs[,j][cov2$treatment=='2'],col="red", pch=1)
# Removing all SVs still results in the PC1 strcuture
# are the PC1s the same?
cor.test(ADfiftypcs[,1],ADfiftypcs5[,1]) 
# 98% Pearson correlation - yes 
# The structure is not created by the SVs...
# but does it exist before correcting for any batch effects?
ADs.raw.Ctrl <- rowADs(voom.countsm[,Ctrl_flies])
ADs.raw.HS <- rowADs(voom.countsm[,HS_flies])
ADs.raw.merge <- cbind(ADs.raw.Ctrl,ADs.raw.HS)
ADs.raw <- ADs.raw.merge[,colnames(voom.countsm)]
ADpca.raw <- prcomp(t(ADs.raw))
plot(ADpca.raw$x[,1] ~ ADpca.raw$x[,2])
ADpca.raw.w = NULL
for (i in 1:length(ADpca.raw$sdev)) {
  ADpca.raw.w[i]= (ADpca.raw$sdev[i])^2/sum(ADpca.raw$sdev^2) 
}
plot(ADpca.raw.w[1:50])
ADfiftypcs.raw <- ADpca.raw$x[, 1:50]
par(mfrow=c(1,2))
i = 2
j = 1
plot(ADfiftypcs.raw[,2]~ ADfiftypcs.raw[,1], col="white", pch=16, main='No batch correction',ylab = "Abs Deviation PC2",xlab = "Abs Deviation PC1")
points(ADfiftypcs.raw[,i][cov2$treatment=="1"]~ ADfiftypcs.raw[,j][cov2$treatment=='1'], col="blue", pch=1)
points(ADfiftypcs.raw[,i][cov2$treatment=="2"]~ ADfiftypcs.raw[,j][cov2$treatment=='2'],col="red", pch=1)
# Removing 1 SV still results in the PC1 structure
plot(fiftypcs.raw[,2]~ fiftypcs.raw[,1], col="white", pch=16, main="Count-based PCA before correction",ylab = "PC2",xlab = "PC1")
points(fiftypcs.raw[,i][cov2$treatment=="1"]~ fiftypcs.raw[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs.raw[,i][cov2$treatment=="2"]~ fiftypcs.raw[,j][cov2$treatment=='2'],col="red", pch=1)
# Removing all SVs still results in the PC1 strcuture
# are the PC1s the same?
cor.test(ADfiftypcs.raw[,1],ADfiftypcs5[,1]) 
# 97% Pearson correlation - they are the same
# the structure exists before and after correcting for the modelled mean effect of day
# Therefore, not a batch correction artefact

# Is there a count PC at any point of the procedure that correlates with it?
cor.test(ADfiftypcs5[,1],fiftypcs[,1]) # -0.0037 - nothing 
cor.test(ADfiftypcs5[,1],fiftypcs5[,1]) # 20% - nothing
cor.test(ADfiftypcs5[,1],fiftypcs.raw[,1]) 
# 80% Pearson correlation - matches the PC1 before correcting for batches
# What does this mean? 
# That the day effect on the variance is very strong?

# 5) Are there other weird PCs?
for (PCnum in c(1,3,5)){
par(mfcol=c(1,1))
j= PCnum
i= PCnum+1
plot(ADfiftypcs5[,i]~ ADfiftypcs5[,j], col="white", pch=16,ylab = paste0("Abs Deviations PC",i),xlab =  paste0("Abs Deviations PC",j))
points(ADfiftypcs5[,i][cov2$treatment=="1"]~ ADfiftypcs5[,j][cov2$treatment=='1'], col="blue", pch=1)
points(ADfiftypcs5[,i][cov2$treatment=="2"]~ ADfiftypcs5[,j][cov2$treatment=='2'],col="red", pch=1)
}
# They are a small portion of the variance 

##########################################
##### End) Save the Work Environment #####
##########################################

save.image(file='ExplainingDevationsPlot.RData')
