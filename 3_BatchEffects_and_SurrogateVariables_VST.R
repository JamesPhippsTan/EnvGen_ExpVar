# Identifying surrogate variables in dataset across different transformations
# Transcriptome datasets of Adult Netherlands D.mel heads on normal vs high sugar diets
# The goal is to infer hidden technical confounding variables present in each transformed expression dataset

# Script by James Tan
# Adapted from Script 'MakingGeneExpressionMatrix_head_HS&CTRL.R' and 'ExploringBatchEffects' by Dr. Luisa Pallares 
# Last Updated: 6/8/24


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


# Working Directory 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")

# Load saved script environment

#######################
##### Datasets ########
#######################

# Metadata on all samples, e.g., conditions and batches
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

nullmodel <- model.matrix(~treatment)

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
# Run anova after each step to assess if there are any residual effects


###############################
##### (2) Transformation ######
###############################

# VST transform count data
data <- DESeqDataSetFromMatrix(countData = raw.counts,
                               colData = cov2,
                               design =  ~treatment) 
vst.data = vst(data, blind=F) 
vst.countsm <- as.matrix(assay(vst.data))


########################################################
##### (3) Remaining batch effect correction ############
########################################################

# Compare counts before and after correcting for batch effects
# Run PCA and plot PCs to visually confirm that each batch effect has been corrected for 
# If you cannot see unique batch-specific clusters in the PCA plots, then the correction has done its job

#############################
##### (3a) VST counts  ######
#############################

pca.raw <- prcomp(t(vst.countsm))
pca.raw.w = NULL
for (i in 1:length(pca.raw$sdev)) {
  pca.raw.w[i]= (pca.raw$sdev[i])^2/sum(pca.raw$sdev^2) 
}
plot(pca.raw.w[1:50])
fiftypcs.raw <- pca.raw$x[, 1:50]

# Check anova against batch effects together
m.raw <-lm(pca.raw$x[,1:100] ~ treatment + rnalib + rnaseq + egg + plate + well)
anova(m.raw)
#Type II MANOVA Tests: Pillai test statistic
#Df test stat approx F num Df den Df    Pr(>F)    
#treatment  1    0.7632  26.8358    200   1665 < 2.2e-16 ***
#  well      95   10.6622   1.1119  19000 167105 < 2.2e-16 ***
#  rnalib     2    1.1178  10.5555    400   3332 < 2.2e-16 ***
#  rnaseq     4    2.3572  11.9662    800   6672 < 2.2e-16 ***
#  egg        2    1.1629  11.5714    400   3332 < 2.2e-16 ***
#  plate      6    3.7938  14.3589   1200  10020 < 2.2e-16 ***
# Well approximate F is very small, regardless of how many PCs are modelled (I tried 50 - 100)
# Small enough to be ignored.

# Plot the PCs for variables.
i = 3
j = 2

#plot treatment 
par(mfrow=c(1,1))
plot(fiftypcs.raw[,i]~ fiftypcs.raw[,j], col="white", pch=16, main="Treatment")
points(fiftypcs.raw[,i][cov2$treatment=="1"]~ fiftypcs.raw[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs.raw[,i][cov2$treatment=="2"]~ fiftypcs.raw[,j][cov2$treatment=='2'],col="red", pch=1)

# Plot well - Very messy and cannot visually inspect. This was not continued further.
plot(fiftypcs.raw[,i]~ fiftypcs.raw[,j], col="white", pch=16, main="Well")
unique_wells <- unique(well)
c_wells <- sample(colours(),length(unique_wells))
for (well_number in 1:length(unique_wells)){
  points(fiftypcs.raw[,i][cov2$well==unique_wells[well_number]]~ fiftypcs.raw[,j][cov2$well==unique_wells[well_number]], col=c_wells[well_number], pch=1)
}

# Plot sequencing batch
plot(fiftypcs.raw[,i] ~ fiftypcs.raw[,j], col="white", pch=16,  main="Seq batch")
points(fiftypcs.raw[,i][cov2$RNAseqBatch==1]~ fiftypcs.raw[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(fiftypcs.raw[,i][cov2$RNAseqBatch==2]~ fiftypcs.raw[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(fiftypcs.raw[,i][cov2$RNAseqBatch==3]~ fiftypcs.raw[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(fiftypcs.raw[,i][cov2$RNAseqBatch==4]~ fiftypcs.raw[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(fiftypcs.raw[,i][cov2$RNAseqBatch==5]~ fiftypcs.raw[,j][cov2$RNAseqBatch==5],col="grey", pch=1)

#plot library prep batch 
plot(fiftypcs.raw[,i] ~ fiftypcs.raw[,j], col="white", pch=16,  main="Lib prep batch")
points(fiftypcs.raw[,i][cov2$RNAlibBatch==1]~ fiftypcs.raw[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(fiftypcs.raw[,i][cov2$RNAlibBatch==2]~ fiftypcs.raw[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(fiftypcs.raw[,i][cov2$RNAlibBatch==3]~ fiftypcs.raw[,j][cov2$RNAlibBatch==3],col="black", pch=1)

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

#plot plating day batch
plot(fiftypcs.raw[,i] ~ fiftypcs.raw[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[1]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[2]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[3]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[4]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[5]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[6]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
points(fiftypcs.raw[,i][cov2$platingBatch==ee[7]]~ fiftypcs.raw[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)


############################################
##### (3b) Counts - without well only ######
############################################

covariates_nowell <- model.matrix(~well)
counts.adj_nowell <- limma::removeBatchEffect(vst.countsm, covariates = covariates_nowell[,-1],  design = nullmodel) 
pca_nowell <- prcomp(t(counts.adj_nowell))
pca_nowell.w = NULL
for (i in 1:length(pca_nowell$sdev)) {
  pca_nowell.w[i]= (pca_nowell$sdev[i])^2/sum(pca_nowell$sdev^2) 
}
plot(pca_nowell.w[1:50]) # PC 1 still 20% of variance

fiftypcs_nowell <- pca_nowell$x[, 1:50]
par(mfrow=c(1,1))

# Let's look at the 1st and 2nd PCs
i = 2
j = 1

#plot treatment - this effect should NOT be lost
plot(fiftypcs_nowell[,i]~ fiftypcs_nowell[,j], col="white", pch=16, main="Treatment")
points(fiftypcs_nowell[,i][cov2$treatment=="1"]~ fiftypcs_nowell[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs_nowell[,i][cov2$treatment=="2"]~ fiftypcs_nowell[,j][cov2$treatment=='2'],col="red", pch=1)

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

# Well hardly changes the first two PCs of the distribution


##################################################
##### (3c) Counts - batch effects with well ######
##################################################

covariates <- model.matrix(~rnalib+rnaseq+egg+plate)
counts.adj <- limma::removeBatchEffect(vst.countsm, covariates = covariates[,-1],  design = nullmodel) 
pca <- prcomp(t(counts.adj))
pca.w = NULL
for (i in 1:length(pca$sdev)) {
  pca.w[i]= (pca$sdev[i])^2/sum(pca$sdev^2) 
}
plot(pca.w[1:50]) # PC1 explains 12% of variance

fiftypcs <- pca$x[, 1:50]
par(mfrow=c(2,2))

# Check anova against treatment and batch effects together
m <-lm(pca$x[,1:200] ~ treatment + well + rnalib + rnaseq + egg + plate)
Anova(m) # effects have been regressed out

# Let's look at the 1st and 2nd PCs
i = 2
j = 1

# plot treatment - this effect should NOT be lost
plot(fiftypcs[,i]~ fiftypcs[,j], col="white", pch=16, main="Treatment")
points(fiftypcs[,i][cov2$treatment=="1"]~ fiftypcs[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs[,i][cov2$treatment=="2"]~ fiftypcs[,j][cov2$treatment=='2'],col="red", pch=1)

# plot sequencing batch - is even now, but the data is still structured at the first pc
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

# plot egglay batch - is even now, but the data is still structured at the first pc
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


##########################################################
##### (3d) Counts - batch effects and well removed #######
##########################################################

covariates_withwell <- model.matrix(~rnalib+rnaseq+egg+plate+well)
counts.adj_withwell <- limma::removeBatchEffect(vst.countsm, covariates = covariates_withwell[,-1],  design = nullmodel) 
pca_withwell <- prcomp(t(counts.adj_withwell))
pca_withwell.w = NULL
for (i in 1:length(pca_withwell$sdev)) {
  pca_withwell.w[i]= (pca_withwell$sdev[i])^2/sum(pca_withwell$sdev^2) 
}
plot(pca_withwell.w[1:50]) 

fiftypcs_withwell <- pca_withwell$x[, 1:50]
par(mfrow=c(1,1))

# Let's look at the 1st and 2nd PCs
i = 2
j = 1

#plot treatment - this effect should NOT be lost
plot(fiftypcs_withwell[,i]~ fiftypcs_withwell[,j], col="white", pch=16, main="Treatment")
points(fiftypcs_withwell[,i][cov2$treatment=="1"]~ fiftypcs_withwell[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs_withwell[,i][cov2$treatment=="2"]~ fiftypcs_withwell[,j][cov2$treatment=='2'],col="red", pch=1)
# There is quite a bit of PC1 structure still present
# Less than if well isn't included

#plot sequencing batch - is even now, but the data is still structured at the first pc
plot(fiftypcs_withwell[,i] ~ fiftypcs_withwell[,j], col="white", pch=16,  main="Seq batch")
points(fiftypcs_withwell[,i][cov2$RNAseqBatch==1]~ fiftypcs_withwell[,j][cov2$RNAseqBatch==1], col="blue", pch=1)
points(fiftypcs_withwell[,i][cov2$RNAseqBatch==2]~ fiftypcs_withwell[,j][cov2$RNAseqBatch==2],col="red", pch=1)
points(fiftypcs_withwell[,i][cov2$RNAseqBatch==3]~ fiftypcs_withwell[,j][cov2$RNAseqBatch==3],col="green", pch=1)
points(fiftypcs_withwell[,i][cov2$RNAseqBatch==4]~ fiftypcs_withwell[,j][cov2$RNAseqBatch==4],col="pink", pch=1)
points(fiftypcs_withwell[,i][cov2$RNAseqBatch==5]~ fiftypcs_withwell[,j][cov2$RNAseqBatch==5],col="grey", pch=1)

#plot library prep batch - is even now, but the data is still structured at the first pc
plot(fiftypcs_withwell[,i] ~ fiftypcs_withwell[,j], col="white", pch=16,  main="Lib prep batch")
points(fiftypcs_withwell[,i][cov2$RNAlibBatch==1]~ fiftypcs_withwell[,j][cov2$RNAlibBatch==1], col="blue", pch=1)
points(fiftypcs_withwell[,i][cov2$RNAlibBatch==2]~ fiftypcs_withwell[,j][cov2$RNAlibBatch==2],col="red", pch=1)
points(fiftypcs_withwell[,i][cov2$RNAlibBatch==3]~ fiftypcs_withwell[,j][cov2$RNAlibBatch==3],col="black", pch=1)

#plot egglay batch - is even now, but the data is still structured at the first pc
plot(fiftypcs_withwell[,i] ~ fiftypcs_withwell[,j], col="white", pch=16,  main="Egglay batch")
c <- sample(colours(),3)
for (e in 1){
  points(fiftypcs_withwell[,i][cov2$egglayBatch==e]~ fiftypcs_withwell[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 2){
  points(fiftypcs_withwell[,i][cov2$egglayBatch==e]~ fiftypcs_withwell[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}
for (e in 3){
  points(fiftypcs_withwell[,i][cov2$egglayBatch==e]~ fiftypcs_withwell[,j][cov2$egglayBatch==e], col=c[e], pch=16)
}

#plot plating day batch - is even now, but the data is still structured at the first pc
plot(fiftypcs_withwell[,i] ~ fiftypcs_withwell[,j], col="white", pch=16,  main="Plating day batch")
c <- sample(colours(),7)
ee <- c(21:27)
points(fiftypcs_withwell[,i][cov2$platingBatch==ee[1]]~ fiftypcs_withwell[,j][cov2$platingBatch==ee[1]], col=c[1], pch=1)
points(fiftypcs_withwell[,i][cov2$platingBatch==ee[2]]~ fiftypcs_withwell[,j][cov2$platingBatch==ee[2]], col=c[2], pch=1)
points(fiftypcs_withwell[,i][cov2$platingBatch==ee[3]]~ fiftypcs_withwell[,j][cov2$platingBatch==ee[3]], col=c[3], pch=1)
points(fiftypcs_withwell[,i][cov2$platingBatch==ee[4]]~ fiftypcs_withwell[,j][cov2$platingBatch==ee[4]], col=c[4], pch=1)
points(fiftypcs_withwell[,i][cov2$platingBatch==ee[5]]~ fiftypcs_withwell[,j][cov2$platingBatch==ee[5]], col=c[5], pch=1)
points(fiftypcs_withwell[,i][cov2$platingBatch==ee[6]]~ fiftypcs_withwell[,j][cov2$platingBatch==ee[6]], col=c[6], pch=1)
points(fiftypcs_withwell[,i][cov2$platingBatch==ee[7]]~ fiftypcs_withwell[,j][cov2$platingBatch==ee[7]], col=c[7], pch=1)
# Well hardly changes the first two PCs of the distribution

# There is remaining structure in the data, use sva to find surrogate batch effects


##################################
##### (4) SV identification ######
##################################

# Finding surrogate variables - tried a range from 1 - 10
mod <-  model.matrix(~ rnalib+rnaseq+egg+plate+well+treatment)
mod0 <- model.matrix(~ rnalib+rnaseq+egg+plate+well)
sva.sva4 <- sva(vst.countsm, mod, mod0,n.sv=10)
sv1.4<- as.numeric(sva.sva4$sv[,1])

# Collecting other svs to see if they may reduce hidden structure in the data further
sv2.4 <- as.numeric(sva.sva4$sv[,2])
sv3.4 <- as.numeric(sva.sva4$sv[,3])
sv4.4 <- as.numeric(sva.sva4$sv[,4])
sv5.4 <- as.numeric(sva.sva4$sv[,5])
sv6.4 <- as.numeric(sva.sva4$sv[,6])


##############################################
##### (3b) Counts - batch effects - 1 SV #####
##############################################

covariates2 <- model.matrix(~rnalib+rnaseq+egg+plate+well+sv1.4)
counts.adj2 <- limma::removeBatchEffect(vst.countsm, covariates = covariates2[,-1],  design = nullmodel) 
pca2 <- prcomp(t(counts.adj2))
plot(pca2$x[,1] ~ pca2$x[,2])
pca2.w = NULL
for (i in 1:length(pca2$sdev)) {
  pca2.w[i]= (pca2$sdev[i])^2/sum(pca2$sdev^2) 
}
plot(pca2.w[1:50])
fiftypcs2 <- pca2$x[, 1:50]

# Check for residual batch effects
mcounts.adj2 <-lm(pca2$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate +well+ sv1.4 + sv2.4 + sv2.4 + sv4.4 + sv5.4)
mcounts.adj2_well <-lm(pca2$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate+ sv1.4 + sv2.4 + sv2.4 + sv4.4 + sv5.4)
anova(mcounts.adj2)
anova(mcounts.adj2_well)
anova(mcounts.adj2,mcounts.adj2_well)
par(mfrow=c(1,1))

# Let's look at the 1st and 2nd PCs
i = 2
j = 1

#plot treatment - this effect should NOT be lost
plot(fiftypcs2[,i]~ fiftypcs2[,j], col="white", pch=16, main="Treatment")
points(fiftypcs2[,i][cov2$treatment=="1"]~ fiftypcs2[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs2[,i][cov2$treatment=="2"]~ fiftypcs2[,j][cov2$treatment=='2'],col="red", pch=1)
# Large PC1 and PC2 structure still present

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


##############################################
##### (3c) Counts - batch effects - 2 SVs ####
##############################################

# Correction for batch effects and 2 surrogate variables
covariates3 <- model.matrix(~rnalib+rnaseq+egg+plate+well+sv1.4+sv2.4)
counts.adj3 <- limma::removeBatchEffect(vst.countsm, covariates = covariates3[,-1],  design = nullmodel) 

pca3 <- prcomp(t(counts.adj3))
plot(pca3$x[,1] ~ pca3$x[,2])
pca3.w = NULL
for (i in 1:length(pca3$sdev)) {
  pca3.w[i]= (pca3$sdev[i])^2/sum(pca3$sdev^2) 
}
plot(pca3.w[1:50])
fiftypcs3 <- pca3$x[, 1:50]

# Check for residual batch effects
mcounts.adj3 <-lm(pca3$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate+ sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4)
mcounts.adj3_well <-lm(pca3$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate +well+ sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4)
anova(mcounts.adj3)
anova(mcounts.adj3_well)
anova(mcounts.adj3,mcounts.adj3_well)

par(mfrow=c(1,1))

# Let's look at the 1st and 2nd PCs
j = 1
i = 2

#plot treatment - this effect should NOT be lost
plot(fiftypcs3[,i]~ fiftypcs3[,j], col="white", pch=16, main="Treatment")
points(fiftypcs3[,i][cov2$treatment=="1"]~ fiftypcs3[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs3[,i][cov2$treatment=="2"]~ fiftypcs3[,j][cov2$treatment=='2'],col="red", pch=1)
# PC2 structure stil present, PC1 structure greatly reduced

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


##############################################
##### (3d) Counts - batch effects - 3 SVs ####
##############################################

# Correction for batch effects and 3 surrogate variables
covariates4 <- model.matrix(~rnalib+rnaseq+egg+plate+well+sv1.4+sv2.4+sv3.4)
counts.adj4 <- limma::removeBatchEffect(vst.countsm, covariates = covariates4[,-1],  design = nullmodel) 

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

#plot treatment - this effect should NOT be lost
plot(fiftypcs4[,i]~ fiftypcs4[,j], col="white", pch=16, main="Treatment")
points(fiftypcs4[,i][cov2$treatment=="1"]~ fiftypcs4[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs4[,i][cov2$treatment=="2"]~ fiftypcs4[,j][cov2$treatment=='2'],col="red", pch=1)
# PC1 and PC2 structure greatly diminished
# Create some mild PC3 structure

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


##############################################
##### (3e) Counts - batch effects - 4 SVs ####
##############################################

# Correction for batch effects and 4 surrogate variables
covariates5 <- model.matrix(~rnalib+rnaseq+egg+plate+well+sv1.4+sv2.4+sv3.4+sv4.4)
counts.adj5 <- limma::removeBatchEffect(vst.countsm, covariates = covariates5[,-1],  design = nullmodel) 

pca5 <- prcomp(t(counts.adj5))
plot(pca5$x[,1] ~ pca5$x[,2])
pca5.w = NULL
for (i in 1:length(pca5$sdev)) {
  pca5.w[i]= (pca5$sdev[i])^2/sum(pca5$sdev^2) 
}
plot(pca5.w[1:50])
fiftypcs5 <- pca5$x[, 1:50]

# Check for residual batch effects
mcounts.adj5 <-lm(pca5$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4)
mcounts.adj5_well <-lm(pca5$x[,1:200] ~ treatment + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4 + well)
anova(mcounts.adj5)
anova(mcounts.adj5_well)
anova(mcounts.adj5,mcounts.adj5_well)

par(mfrow=c(1,1))


# Let's look at the 1st and 2nd PCs
i = 3
j = 2

#plot treatment - this effect should NOT be lost
plot(fiftypcs5[,i]~ fiftypcs5[,j], col="white", pch=16, main="Treatment")
points(fiftypcs5[,i][cov2$treatment=="1"]~ fiftypcs5[,j][cov2$treatment=='1'], col="blue", pch=1)
points(fiftypcs5[,i][cov2$treatment=="2"]~ fiftypcs5[,j][cov2$treatment=='2'],col="red", pch=1)
# Creates some new PC2 structure (points beyond the 0)
# PC3 structure introduced

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
# No further correction needed


#######################################
##### Checking effect of SVs ##########
#######################################

# With well included
mw1 <-lm(pca.raw$x[,1:200] ~ treatment + well + rnalib + rnaseq + egg + plate)
mw2 <-lm(pca.raw$x[,1:200] ~ treatment + well + rnalib + rnaseq + egg + plate + sv1.4)
mw3 <-lm(pca.raw$x[,1:200] ~ treatment + well + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4)
mw4 <-lm(pca.raw$x[,1:200] ~ treatment + well + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4)
mw5 <-lm(pca.raw$x[,1:200] ~ treatment + well + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4)
mw6 <-lm(pca.raw$x[,1:200] ~ treatment + well + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4)
mw7 <-lm(pca.raw$x[,1:200] ~ treatment + well + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4 + sv6.4)
anova(mw1,mw2,mw3,mw4,mw5,mw6,mw7)

# Without well 
m1 <-lm(pca.raw$x[,1:200] ~ treatment  + rnalib + rnaseq + egg + plate)
m2 <-lm(pca.raw$x[,1:200] ~ treatment  + rnalib + rnaseq + egg + plate + sv1.4)
m3 <-lm(pca.raw$x[,1:200] ~ treatment  + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4)
m4 <-lm(pca.raw$x[,1:200] ~ treatment  + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4)
m5 <-lm(pca.raw$x[,1:200] ~ treatment  + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4)
m6 <-lm(pca.raw$x[,1:200] ~ treatment  + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4)
m7 <-lm(pca.raw$x[,1:200] ~ treatment  + rnalib + rnaseq + egg + plate + sv1.4 + sv2.4 + sv3.4 + sv4.4 + sv5.4 + sv6.4)
anova(m1,m2,m3,m4,m5,m6,m7)

# Compare all models
anova(m1,mw1)
anova(m1,m2,m3,m4,m5,m6,m7,mw1,mw2,mw3,mw4,mw5,mw6,mw7)


# The well models are always significantly better
# Adding more surrogate variables increases the fit of the model
# However, the PC1 and PC2 structure issue is solved within the first three PCs
# Therefore, I will proceed with the 3 PCs


##########################################
##### End) Save the Work Environment #####
##########################################

write.table(sv1.4, 'VST_sv1_10.txt')
write.table(sv2.4, 'VST_sv2_10.txt')
write.table(sv3.4, 'VST_sv3_10.txt')
write.table(sv4.4, 'VST_sv4_10.txt')
write.table(sv5.4, 'VST_sv5_10.txt')

save.image(file='NvsHS_BatchEffects_and_SurrogateVariables_VST.RData')
