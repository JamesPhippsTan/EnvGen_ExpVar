# Transcriptome variability analyses of Adult Netherlands D.mel heads on normal vs high sugar diets
# Tests for differences in mean and variabilities between conditions

# Script by James Tan

# Includes other methods for compariosn
# BF test for different variabilities with VST for VE analysis

# BF's test section adapted from (https://rpubs.com/LiYumei/806213)

# Last Updated: 15/9/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(DESeq2)
library(dplyr)
library(edgeR)
library(limma)
library(tidyr)
library(tibble)
library(car)
library(matrixStats)
library(ggplot2)
library(ggrepel)
library(vGWAS)
library(gamlss)
library(ggvenn)

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV")
load(file='DE_DV_BF_Test.RData') 

# Obtain functions used to calculate adjusted variabilities
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('Variability_Functions.R')


#####################
##### Datasets ######
#####################

# Metadata on all samples, e.g., conditions and batches
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)

# Surrogate variables, estimated in the NvsHS_SurrogateVariables R scripts
# Based on manual PC inspection - Include 3 SVs for VST; 4 SVs for voom
tmm.voom.sv1 <- read.table('TMM_Voom_sv1_10.txt')
tmm.voom.sv1 <- tmm.voom.sv1$x
tmm.voom.sv2 <- read.table('TMM_Voom_sv2_10.txt')
tmm.voom.sv2 <- tmm.voom.sv2$x
tmm.voom.sv3 <- read.table('TMM_Voom_sv3_10.txt')
tmm.voom.sv3 <- tmm.voom.sv3$x
tmm.voom.sv4 <- read.table('TMM_Voom_sv4_10.txt')
tmm.voom.sv4 <- tmm.voom.sv4$x

# Expression - raw counts of filtered samples (see MakingGeneExpressionMatrix_head_HS&CTRL.R for filtering criteria)
raw.counts <- read.table("RawCounts_noY_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt",h=T,check.names = F)

# Eliminate bimodal genes
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
Bimodal_genes <- read.csv('Bimodal_Genes.csv',row.names = 1) 
Non_bimodal_genes <- setdiff(rownames(raw.counts),rownames(Bimodal_genes))
raw.counts <- raw.counts[Non_bimodal_genes,]

# Only metadata on samples with counts
cov <- tibble::column_to_rownames(head.info, "id")
cov2 <- cov[c(colnames(raw.counts)),]

# Merge expression and metadata 
raw.counts.t <- data.frame(t(raw.counts)) 
raw.counts.t <- tibble::rownames_to_column(raw.counts.t, "id")
head.data <- merge(raw.counts.t,head.info,by = 'id')
head.data <- tibble::column_to_rownames(head.data, "id")


############################
##### Transformations ######
############################

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


##########################
##### TMM Voom ###########
##########################

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


####################################################################################
##### 1) Is the variability significantly different between conditions? ############
####################################################################################

######################################
######### BF test on Voom data #######
######################################

# Perhaps closer to the GAMLSS data, given that the data isn't corrected for the transcriptome-wide mean-variability trend

pvalues.voom <- sapply(1:nrow(voom.counts.bc), function(i){
  data <- cbind.data.frame(gene = as.numeric(t(voom.counts.bc[i,])), conditions)
  result <- leveneTest(gene~conditions, data, center='median')
  p <- result$`Pr(>F)`[1]
  return(p)
})
bh.voom <- p.adjust(pvalues.voom, method = "BH")

# Obtain voom count-based means and MADs for downstream correlations
CtrlMeanVoom <- rowMeans(voomdataN)
HSMeanVoom <- rowMeans(voomdataHS)
MAD_N_Voom <- matrixStats::rowMads(as.matrix(voomdataN))/1.4826
MAD_HS_Voom <- matrixStats::rowMads(as.matrix(voomdataHS))/1.4826

# Output test results 
BFTestVoomResult <- data.frame(pval.BFtest.Voom = pvalues.voom,
                               padj.BFtest.Voom = bh.voom,
                               Ctrl.Mean.Voom = CtrlMeanVoom, 
                               HS.Mean.Voom = HSMeanVoom,
                               Ctrl.MAD.Voom = MAD_N_Voom, 
                               HS.MAD.Voom = MAD_HS_Voom,
                               logFC.MAD = MAD_HS_Voom/MAD_N_Voom)
rownames(BFTestVoomResult) <- rownames(voom.counts.bc)
BFTestVoomResult$Gene_id <- rownames(voom.counts.bc)

# How many significant results?
BF_sig <- subset(BFTestVoomResult,padj.BFtest.Voom<0.05) # 7231
BF_inc <- subset(BFTestVoomResult,padj.BFtest.Voom<0.05&logFC.MAD>1) # 6868
BF_dec <- subset(BFTestVoomResult,padj.BFtest.Voom<0.05&logFC.MAD<1) # 363

# Compare with GAMLSS
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_GeneSets/")
GAMLSS_sig <- read.csv('ADVGs.csv',row.names = 1) # 7685
GAMLSS_inc <- subset(GAMLSS_sig,logFC.BCV>0) # 6896
GAMLSS_dec <- subset(GAMLSS_sig,logFC.BCV<0) # 789

# Intersections
# With the logCPM data
length(intersect(rownames(BFvoom_sig),rownames(GAMLSS_sig))) # 7030 -> 91% overlap with GAMLSS
length(intersect(rownames(BFvoom_inc),rownames(GAMLSS_inc))) # 6312 -> 91% overlap with GAMLSS
length(intersect(rownames(BFvoom_dec),rownames(GAMLSS_dec))) # 515 -> 65% overlap with GAMLSS
# 2 papers use this method...(unclear for other papers whether log was used)

# Merge with gamlss dataframe and save
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\")
gamlss_NB_final <- read.csv('AllGenes_DE_DV.csv',row.names = 1) # 7685

# Only need 1st 16 columns
DV_testresults <- merge(gamlss_NB_final[,1:16],BFTestVoomResult,all=T,by='Gene_id')
View(DV_testresults)
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\")
write.csv(DV_testresults,'Both_DV_test_results.csv') # 7685


##########################################################
##### Variability change by mean plots ###################
##########################################################

# Get DVG overlap
gene_set1 <- rownames(GAMLSS_sig)
gene_set2 <- rownames(BFvoom_sig)
gene_lists <- list(
  "GAMLSS" = rownames(GAMLSS_sig),
  "Brown-Forsythe test" = rownames(BFvoom_sig)
)

p <- ggvenn(gene_lists,
            fill_color = c("#66c2a5", "#fc8d62"),
            stroke_size = 0.4,
            set_name_size = 4,
            text_size = 4)
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\")
p
ggsave(plot = p,filename = 'GAMLSS_BF_Overlap.svg',width=4,height=4,dpi=300)

# Estimated mean by estimated BCV
BCV_Mean_Ctrl_plot <- ggplot(data=gamlss_NB_final, aes(x=rank(cpm.ctrl), y=BCV.ctrl)) + 
  geom_point(alpha = 0.5,col='darkgrey') +
  xlab("Ctrl mean estimated CPM rank") + 
  ylab("Ctrl estimated BCV") + theme_classic()+
  theme(legend.position = "none")
cor.test(gamlss_NB_final$cpm.ctrl,gamlss_NB_final$BCV.ctrl,method = 'spearman')
BCV_Mean_Ctrl_plot

# Estimated mean by estimated BCV
BCV_Mean_Ctrl_plot <- ggplot(data=gamlss_NB_final, aes(x=rank(cpm.ctrl), y=BCV.ctrl)) + 
  geom_point(alpha = 0.5,col='darkgrey') +
  xlab("Ctrl mean estimated CPM rank") + 
  ylab("Ctrl estimated BCV") + theme_classic()+
  theme(legend.position = "none")
cor.test(gamlss_NB_final$cpm.ctrl,gamlss_NB_final$BCV.ctrl,method = 'spearman')
BCV_Mean_Ctrl_plot

BCV_Mean_HS_plot <- ggplot(data=gamlss_NB_final, aes(x=rank(cpm.hs), y=BCV.hs)) + 
  geom_point(alpha = 0.5,col='darkgrey') +
  xlab("HS mean estimated CPM rank") + 
  ylab("HS estimated BCV") + theme_classic()+
  theme(legend.position = "none")
cor.test(gamlss_NB_final$cpm.hs,gamlss_NB_final$BCV.hs,method = 'spearman')
BCV_Mean_HS_plot

# Variability change by mean rank
DelBCV_Mean_HS_plot <- ggplot(data=gamlss_NB_final, aes(x=rank(cpm.hs), y=log2FC.BCV, col=Different_Variability)) + 
  geom_point(alpha = 0.5) +
  ylab(expression(log[2]~"FC in variability (HS BCV / Ctrl BCV)")) + 
  xlab("HS mean estimated CPM rank") + theme_classic()+
  scale_color_manual(
    values = c("No" = "darkgrey", "Yes" = "#611BB8"))+
  theme(legend.position = "none")
cor.test(gamlss_NB_final$cpm.hs,gamlss_NB_final$log2FC.BCV,method = 'spearman')
# Small positive correlation
DelBCV_Mean_HS_plot

# Variability change by mean change
DelBCV_DelMean_Ctrl_plot <- ggplot(data=gamlss_NB_final, aes(x=log2FC.cpm, y=log2FC.BCV, col=Different_Variability)) + 
  geom_point(alpha = 0.5) +
  ylab(expression(log[2]~"FC in variability (HS BCV / Ctrl BCV)")) + 
  xlab(expression(log[2]~"FC in mean CPM (HS / Ctrl)")) + theme_classic()+
  scale_color_manual(
    values = c("No" = "darkgrey", "Yes" = "#611BB8"))+
  theme(legend.position = "none")
cor.test(gamlss_NB_final$log2FC.cpm,gamlss_NB_final$log2FC.BCV,method = 'spearman')
# Small negative correlation
DelBCV_DelMean_Ctrl_plot

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\")
ggsave(plot = BCV_Mean_Ctrl_plot,filename = 'BCV_Mean_Ctrl_plot.svg',width=3.5,height=3.5,dpi=300)
ggsave(plot = BCV_Mean_HS_plot,filename = 'BCV_Mean_HS_plot.svg',width=3.5,height=3.5,dpi=300)
ggsave(plot = DelBCV_Mean_HS_plot,filename = 'DelBCV_Mean_HS_plot.svg',width=3.5,height=3.5,dpi=300)
ggsave(plot = DelBCV_DelMean_Ctrl_plot,filename = 'DelBCV_DelMean_Ctrl_plot.svg',width=3.5,height=3.5,dpi=300)

# Voom versions

# Estimated mean by estimated VoomVar
VoomVar_Mean_Ctrl_plot <- ggplot(data=BFTestVoomResult, aes(x=rank(CtrlMeanVoom), y=CtrlVariabilityVoom)) + 
  geom_point(alpha = 0.5,col='darkgrey') +
  xlab("Ctrl mean logCPM rank") + 
  ylab("Ctrl calculated MAD") + theme_classic()+
  theme(legend.position = "none")
cor.test(BFTestVoomResult$CtrlMeanVoom,BFTestVoomResult$CtrlVariabilityVoom,method = 'spearman')
VoomVar_Mean_Ctrl_plot

# Estimated mean by estimated VoomVar
VoomVar_Mean_HS_plot <- ggplot(data=BFTestVoomResult, aes(x=rank(HSMeanVoom), y=HSVariabilityVoom)) + 
  geom_point(alpha = 0.5,col='darkgrey') +
  xlab("HS mean logCPM rank") + 
  ylab("HS calculated MAD") + theme_classic()+
  theme(legend.position = "none")
cor.test(BFTestVoomResult$HSMeanVoom,BFTestVoomResult$HSVariabilityVoom,method = 'spearman')
VoomVar_Mean_HS_plot


# Variability change by mean rank
BFTestVoomResult$Different_Variability_BF <- ifelse(BFTestVoomResult$Var.BH.pval.voom<0.05,'Yes','No')
DelVoomVar_Mean_HS_plot <- ggplot(data=BFTestVoomResult, aes(x=rank(HSMeanVoom), y=log2(Ratio_Var2_Var1_Voom), col=Different_Variability_BF)) + 
  geom_point(alpha = 0.5) +
  ylab(expression(log[2]~"FC in variability (HS MAD / Ctrl MAD)")) + 
  xlab("HS mean logCPM rank") + theme_classic()+
  scale_color_manual(
    values = c("No" = "darkgrey", "Yes" = "#611BB8"))+
  theme(legend.position = "none")
cor.test(BFTestVoomResult$HSMeanVoom,BFTestVoomResult$Ratio_Var2_Var1_Voom,method = 'spearman')
# Small positive correlation
DelVoomVar_Mean_HS_plot

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\")
ggsave(plot = VoomVar_Mean_Ctrl_plot,filename = 'VoomVar_Mean_Ctrl_plot.svg',width=3.5,height=3.5,dpi=300)
ggsave(plot = VoomVar_Mean_HS_plot,filename = 'VoomVar_Mean_HS_plot.svg',width=3.5,height=3.5,dpi=300)
ggsave(plot = DelVoomVar_Mean_HS_plot,filename = 'DelVoomVar_Mean_HS_plot.svg',width=3.5,height=3.5,dpi=300)



# Save workspace
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV")
save.image(file='DE_DV_BF_Test.RData') 


