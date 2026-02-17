# Transcriptome variability analyses of Adult Netherlands D.mel heads on normal vs high sugar diets
# Tests for differences in mean (DE) and variabilities (DV) between conditions

# Script by James Tan

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
library(eulerr)
library(svglite)


# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV")
#load(file='DE_DV_Main.RData')

# Obtain functions for mean-variance plots
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('Variability_Functions.R')

#####################
##### Datasets ######
#####################

# Metadata on all samples, e.g., conditions and batches
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)

# Surrogate variables, estimated in the investigating batch effect R scripts
# Based on manual PC inspection - 4 SVs for voom
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

# Bimodal genes
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
Bimodal_genes <- read.csv('Bimodal_Genes.csv',row.names = 1) 

# Genewise means and variabilities
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
MeanVars <- read.csv('Gene_MeanVar_Table.csv') 

##############################
##### Reformat datasets ######
##############################

# Eliminate bimodal genes
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

# Sums to zero contrast for the batch effects
contrasts(rnalib) <- contr.sum(levels(rnalib)) 
contrasts(rnaseq) <- contr.sum(levels(rnaseq))
contrasts(egg) <- contr.sum(levels(egg))
contrasts(plate) <- contr.sum(levels(plate))
contrasts(well) <- contr.sum(levels(well))

##################################
##### Run GAMLSS DE and DV #######
##################################

raw.counts.DGElist <- DGEList(counts=raw.counts)

# Calculate TMM library size normalisation factors 
raw.counts.DGElist <- calcNormFactors(raw.counts.DGElist, method="TMM")

# Extract TMM-normalised library sizes to explicitly account for this in the model
raw.counts.DGElist$samples$offset <- log(raw.counts.DGElist$samples$lib.size * raw.counts.DGElist$samples$norm.factors)
offset_values <- raw.counts.DGElist$samples$offset

# Full design matrix for model
design <- model.matrix(~rnalib+rnaseq+egg+plate+well+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4+treatment)
# Note: intercept represents the mean of control flies

# GAMLSS function
gamlss_function <- function(gene,sigma.start.estimate,max_cycles) {
  y=raw.counts.DGElist$counts[gene,]
  data <- as.data.frame(cbind(y, design[,-1]))
  data$treatment2 <- as.factor(data$treatment2)
  data$treatment2 <- relevel(data$treatment2, ref = c("0")) # turn into 0s and 1s
  
  # Full_model
  full_model <- as.formula(paste("y ~ ", paste(colnames(data)[-1], collapse = " + "), "+ offset(offset_values)"))
  m0 <- tryCatch(
    gamlss(fo = full_model, sigma.fo = ~treatment, data=data,
           family = NBI(), sigma.start = sigma.start.estimate, n.cyc = max_cycles),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit reduced model by omitting treatment from the estimation of overdispersion: sigma.fo = ~ 1. 
  # Besides not adjusting sigma based on the genome-wide mean-variance correlation, this model corresponds to the GLM model implemented in edgeR.
  m1 <- tryCatch(
    gamlss(fo = full_model, sigma.fo = ~ 1, data=data,
           family = NBI(), sigma.start = sigma.start.estimate, n.cyc = max_cycles),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit reduced model by omitting treatment from the estimation of mean: sigma.fo = ~ 1. 
  notreatment_model <- as.formula(paste("y ~ ", paste(colnames(data)[-c(1, which(colnames(data) == 'treatment2'))], collapse = " + "), "+ offset(offset_values)"))
  m2 <- tryCatch(
    gamlss(fo = notreatment_model, sigma.fo =  ~treatment, data=data,
           family = NBI(), sigma.start = sigma.start.estimate, n.cyc = max_cycles),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Create data frame res to store the results.
  res <- data.frame(
    cpm.ctrl = NA,
    cpm.hs = NA,
    mu.ctrl=NA,
    mu.hs=NA,
    logFC.cpm = NA,
    pval.cpm = NA,
    BCV.ctrl = NA,
    BCV.hs = NA,
    sigma.ctrl=NA,
    sigma.hs=NA,
    logFC.BCV = NA,
    pval.BCV = NA
  )
  
  # Because fitting of the NB model may fail for some genes, check whether all models were fitted successfully. 
  if(!any(sapply(list(m0,m1,m2), is.null))) 
  {
    # Write GAMLSS estimations of gene’s mean (CPM) counts from the m0 model.
    res$cpm.ctrl = exp(m0$mu.coefficients[[1]]+log(1e06))
    res$cpm.hs = exp(m0$mu.coefficients[[1]]+m0$mu.coefficients[["treatment21"]]+log(1e06))
    res$mu.ctrl = m0$mu.coefficients[[1]]
    res$mu.hs = m0$mu.coefficients[[1]]+m0$mu.coefficients[["treatment21"]]
    
    # Calculate log ratio for changes in CPMs between Ctrl and HS flies
    res$logFC.cpm = log2(res$cpm.hs/res$cpm.ctrl)
    
    # GAMLSS log-likelihood ratio (LR) test for a significance of a treatnent effect on gene’s mean (CPM) counts.
    res$pval.cpm = pchisq(2*(logLik(m0)-logLik(m2)), df=m0$df.fit-m2$df.fit, lower=F)[[1]]
    
    # Write GAMLSS estimations of gene’s non-Poisson noise from the m0 model: bcv(μ)=√α.
    res$BCV.ctrl = sqrt(exp(m0$sigma.coefficients[[1]]))
    res$BCV.hs = sqrt(exp(m0$sigma.coefficients[[1]]+m0$sigma.coefficients[["treatment2"]]))
    res$sigma.ctrl = m0$sigma.coefficients[[1]]
    res$sigma.hs = m0$sigma.coefficients[[1]]+m0$sigma.coefficients[["treatment2"]]
    
    # Calculate log ratio for changes in cv(μ) between Ctrl and HS flies
    res$logFC.BCV = log2(res$BCV.hs/res$BCV.ctrl)
    
    # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on non-Poisson noise.
    # padj.BCV – p value of LR test statistic: D_α=-2log⁡[L(μ_j,α_0  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m1 models.
    res$pval.BCV = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
  }
  print(paste0("Gene ",gene," Complete"))
  res
}

# GAMLSS Modelling
gene_id <- 1:nrow(raw.counts.DGElist)
# The default settings are 0.1 as a starting alpha
# As well as 100 cycles maximum to find the best-fitting alpha
gamlss_NB <- lapply(gene_id, gamlss_function, sigma.start.estimate = 0.1, max_cycles=100)
gamlss_NB <- do.call(rbind, gamlss_NB)
rownames(gamlss_NB) <- rownames(raw.counts.DGElist$counts)[gene_id]

# 11 genes had NAs and one has a weird logFC.cpm - p.value profile
gamlss_NA_genes <- rownames(gamlss_NB)[is.na(gamlss_NB[["logFC.cpm"]])]
weird_genes <- c(gamlss_NA_genes,"FBgn0034871")
weird_genes
# See code at the end for an explanation why these genes are weird...
# Repeat the gamlss but with a larger starting parameter of sigma
gamlss_NB_weird_genes <- lapply(weird_genes, gamlss_function, sigma.start.estimate = 1, max_cycles=100)
gamlss_NB_weird_genes <- do.call(rbind, gamlss_NB_weird_genes)
rownames(gamlss_NB_weird_genes) <- weird_genes
View(gamlss_NB_weird_genes)

# Merge the dataframes
gamlss_NB_fine_genes <- setdiff(rownames(gamlss_NB),weird_genes)
gamlss_NB_full <- rbind(gamlss_NB[gamlss_NB_fine_genes,],gamlss_NB_weird_genes)
gamlss_NB_final <- gamlss_NB_full[complete.cases(gamlss_NB_full),]
nrow(gamlss_NB_final) # 4 genes still did not converge

# Apply the BH correction to call significantly increased/decreased genes
gamlss_NB_final$padj.cpm <- p.adjust(gamlss_NB_final$pval.cpm, method = 'BH')
gamlss_NB_final$padj.BCV <- p.adjust(gamlss_NB_final$pval.BCV, method = 'BH')

# Consider raw changes in BCV rather than log changes
gamlss_NB_final$Difference.BCV <- gamlss_NB_final$BCV.hs-gamlss_NB_final$BCV.ctrl
    
# Consider log2 fold changes (default code was logn)
gamlss_NB_final$log2FC.cpm = log2(gamlss_NB_final$cpm.hs/gamlss_NB_final$cpm.ctrl)
gamlss_NB_final$log2FC.BCV = log2(gamlss_NB_final$BCV.hs/gamlss_NB_final$BCV.ctrl)


##############################################
##### 5) Defining categories of DE and VE ####
##############################################

# Save DE genes 
DEGsFinal <- subset(gamlss_NB_final, padj.cpm < 0.05)
DEGsIncreaseHS <- subset(DEGsFinal,logFC.cpm > 0) 
DEGsDecreaseHS <- subset(DEGsFinal,logFC.cpm < 0) 
# Pure DEGs (no variability change)
PDEGsFinal <- subset(gamlss_NB_final, padj.cpm < 0.05 & padj.BCV >= 0.05)
PDEGsIncreaseHS <- subset(PDEGsFinal,logFC.cpm > 0)
nrow(PDEGsIncreaseHS) # 435
PDEGsDecreaseHS <- subset(PDEGsFinal,logFC.cpm < 0)
nrow(PDEGsDecreaseHS) # 96

# Save DV genes
# A(ll)DV genes
ADVGsFinal <- subset(gamlss_NB_final, padj.BCV < 0.05)
ADVGsDecanalizedHS <- subset(ADVGsFinal,logFC.BCV > 0)
nrow(ADVGsDecanalizedHS) #6896
ADVGsCanalizedHS <- subset(ADVGsFinal,logFC.BCV < 0)
nrow(ADVGsCanalizedHS) #789
# Purely DVGs (Call DVG due to convention)
DVGsFinal <- subset(gamlss_NB_final, padj.cpm >= 0.05 & padj.BCV < 0.05)
DVGsDecanalizedHS <- subset(DVGsFinal,logFC.BCV > 0)
nrow(DVGsDecanalizedHS) # 3444
DVGsCanalizedHS <- subset(DVGsFinal,logFC.BCV < 0)
nrow(DVGsCanalizedHS) # 390

# DE and DV genes
DEDVGsFinal <- subset(gamlss_NB_final, padj.cpm < 0.05 & padj.BCV < 0.05)
# DEDVGs that increase vs those that decrease in variance under HS
DEDVGsDecanalizedHS <- subset(DEDVGsFinal,logFC.BCV > 0)
nrow(DEDVGsDecanalizedHS) # 3452
DEDVGsCanalizedHS <- subset(DEDVGsFinal,logFC.BCV < 0)
nrow(DEDVGsCanalizedHS) # 399
# DEDVGs by increased or decreased mean expression
DEDVGsDecanalizedIncreasedHS <- subset(DEDVGsDecanalizedHS,logFC.cpm > 0) 
nrow(DEDVGsDecanalizedIncreasedHS) # 2015 
DEDVGsDecanalizedDecreasedHS <- subset(DEDVGsDecanalizedHS,logFC.cpm < 0) 
nrow(DEDVGsDecanalizedDecreasedHS) # 1437
DEDVGsCanalizedIncreasedHS <- subset(DEDVGsCanalizedHS,logFC.cpm > 0) 
nrow(DEDVGsCanalizedIncreasedHS) # 316
DEDVGsCanalizedDecreasedHS <- subset(DEDVGsCanalizedHS,logFC.cpm < 0) 
nrow(DEDVGsCanalizedDecreasedHS) # 83

# Genes with no change in variability or mean in response to HS
NDE_NDV_Gs <- subset(gamlss_NB_final, padj.cpm >= 0.05 & padj.BCV >= 0.05)
nrow(NDE_NDV_Gs) # 543

#################################################################################
##### (6) Make yes or no labels for all the subcategories of genes  #############
#################################################################################

# Make Gene_id column
gamlss_NB_final$Gene_id <- rownames(gamlss_NB_final)

# All gene sets of interest
gene_sets <- list(Increased_Mean_HS = rownames(DEGsIncreaseHS),
                  Decreased_Mean_HS = rownames(DEGsDecreaseHS),
                 Increased_Variability_HS = rownames(ADVGsDecanalizedHS),
                 Decreased_Variability_HS = rownames(ADVGsCanalizedHS),
                 # Variability by change in mean
                 Increased_Variability_Different_Mean = rownames(DEDVGsDecanalizedHS),
                 Decreased_Variability_Different_Mean = rownames(DEDVGsCanalizedHS),
                 # 9 specific categories last
                 Decreased_Variability_Same_Mean_HS = rownames(DVGsCanalizedHS),
                 Increased_Variability_Same_Mean_HS = rownames(DVGsDecanalizedHS),
                 Decreased_Variability_Decreased_Mean_HS = rownames(DEDVGsCanalizedDecreasedHS),
                 Decreased_Variability_Increased_Mean_HS = rownames(DEDVGsCanalizedIncreasedHS),
                 Increased_Variability_Decreased_Mean_HS = rownames(DEDVGsDecanalizedDecreasedHS),
                 Increased_Variability_Increased_Mean_HS = rownames(DEDVGsDecanalizedIncreasedHS),
                 Same_Variability_Decreased_Mean_HS = rownames(PDEGsDecreaseHS),
                 Same_Variability_Increased_Mean_HS = rownames(PDEGsIncreaseHS),
                 Same_Variability_Same_Mean = rownames(NDE_NDV_Gs)
                )

# Function to return yes or no
check_gene_in_set <- function(gene, gene_set) {
  if (gene %in% gene_set) {
    return("Yes")
  } else {
    return("No")
  }
}

# Apply the function to each gene set to create new columns in the larger dataframe
for (set_name in names(gene_sets)) {
  gene_set <- gene_sets[[set_name]]
  gamlss_NB_final[[set_name]] <- sapply(gamlss_NB_final$Gene_id, check_gene_in_set, gene_set)
}

# Make a table summarizing the number of genes in each category
gene_set_numbers <- data.frame(
  Gene_set = c("Whole Transcriptome",names(gene_sets)),
  Number_of_genes = c(nrow(gamlss_NB_final),sapply(gene_sets, length))
)

# Percentage of transcriptome
format_transcriptome_percent <- function(value){
    round(100*value/8763,digits = 1)
}
gene_set_numbers$Percent_transcriptome <- format_transcriptome_percent(gene_set_numbers$Number_of_genes)

##############################################################
##### (7) Plots using the subcategories of genes  ############
##############################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\")

DV_label_decreased <- paste0(gene_set_numbers['Decreased_Variability_HS','Percent_transcriptome'],"%")
DV_label_increased <- paste0(gene_set_numbers['Increased_Variability_HS','Percent_transcriptome'],"%")
DV_plot <- ggplot(data=gamlss_NB_final, aes(x=log2FC.BCV, y=-log10(pval.BCV), col=Different_Variability)) + 
  geom_point(alpha = 0.5) + geom_hline(yintercept=-log10(max(ADVGsFinal$pval.BCV)), col="gray50",linetype='dashed')+
  xlab(expression(log[2]~"FC in variability (HS BCV / Ctrl BCV)")) + 
  ylab(expression("-"~log[10]~"p-val ("~sigma~")")) + theme_classic()+
  scale_color_manual(
    values = c("No" = "grey", "Yes" = "#611BB8"))+
  theme(legend.position = "none")+ 
  annotate("text", y = 70, x = 2, label = DV_label_increased, size = 5, hjust = 1,vjust = 1)+ 
  annotate("text", y = 70, x = -2, label = DV_label_decreased, size = 5, hjust = 0,vjust = 1)+
  xlim(c(-2,2))
DV_plot
ggsave(plot = DV_plot,filename = "DV_plot.svg",height=3,width=3,dpi=300)

# DE gene volcano plot
DE_label_decreased <- paste0(gene_set_numbers['Decreased_Mean_HS','Percent_transcriptome'],"%")
DE_label_increased <- paste0(gene_set_numbers['Increased_Mean_HS','Percent_transcriptome'],"%")
DE_plot <- ggplot(data=gamlss_NB_final, aes(x=log2FC.cpm, y=-log10(pval.cpm), col=Different_Mean)) + 
  geom_point(alpha = 0.5) + geom_hline(yintercept=-log10(max(DEGsFinal$pval.cpm)), col="gray50",linetype='dashed')+
  xlab(expression(log[2]~"FC in mean CPM (HS / Ctrl)")) + 
  ylab(expression("-"~log[10]~"p-val ("~mu~")")) + theme_classic()+
  scale_color_manual(
    values = c("No" = "grey", "Yes" = "#4F8E4D"))+
  theme(legend.position = "none")+ 
  annotate("text", y = 90, x = 1.5, label = DE_label_increased, size = 5, hjust = 0,vjust = 0.5)+ 
  annotate("text", y = 90, x = -1.5, label = DE_label_decreased, size = 5, hjust = 1,vjust = 0.5)+
  xlim(c(-5,5))
DE_plot
ggsave(plot = DE_plot,filename = "DE_plot.svg",height=3,width=3,dpi=300)


##########################################################
##### (8) Mean plots of the different categories #########
##########################################################

# All unique gene sets against whole transcriptome
gene_sets_to_ridge <- list(All_genes = MeanVars$gene_ID,
                  Same_variability_same_mean = rownames(NDE_NDV_Gs),
                  Lower_variability_lower_mean = rownames(DEDVGsCanalizedDecreasedHS),
                  Lower_variability_higher_mean = rownames(DEDVGsCanalizedIncreasedHS),
                  Higher_variability_lower_mean = rownames(DEDVGsDecanalizedDecreasedHS),
                  Higher_variability_higher_mean = rownames(DEDVGsDecanalizedIncreasedHS),
                  Higher_variability_same_mean = rownames(DVGsDecanalizedHS),
                  Lower_variability_same_mean = rownames(DVGsCanalizedHS),
                  Same_variability_lower_mean = rownames(PDEGsDecreaseHS),
                  Same_variability_higher_mean = rownames(PDEGsIncreaseHS)
)

multiset_gene_ridges <- plot_multiset_gene_ridges(df = MeanVars,
                                                  gene_sets = gene_sets_to_ridge,
                                               gene_column = 'gene_ID',
                                               value_column = 'HS_logCPM_Mean',
                                               x_label = 'HS mean logCPM',
                                               collapse_label = F)
multiset_gene_ridges
ggsave(plot = multiset_gene_ridges,filename = 'All_DEDV_MeanHS_Distributions.svg',width=6,height=6,dpi=300)

# All unique gene sets against whole transcriptome
DEDVGs_to_ridge <- list(All_genes = MeanVars$gene_ID,
                      Higher_mean = rownames(DEGsIncreaseHS),
                      Lower_mean = rownames(DEGsDecreaseHS),
                      Higher_variability = rownames(ADVGsDecanalizedHS),
                      Lower_variability = rownames(ADVGsCanalizedHS)
)
multiset_gene_ridges_DEDVGs <- plot_multiset_gene_ridges(df = MeanVars,
                                                  gene_sets = DEDVGs_to_ridge,
                                                  gene_column = 'gene_ID',
                                                  value_column = 'HS_logCPM_Mean',
                                                  x_label = 'HS mean logCPM',
                                                  collapse_label = T)
multiset_gene_ridges_DEDVGs
ggsave(plot = multiset_gene_ridges_DEDVGs,filename = 'DEDVGs_MeanHS_Distributions.svg',width=3,height=3,dpi=300)


# Interesting gene sets for main figure
gene_sets_to_ridge_subset <- list(Higher_variability = rownames(ADVGsDecanalizedHS),
                                  Lower_variability = rownames(ADVGsCanalizedHS),
                                  All_genes = MeanVars$gene_ID
)

multiset_gene_ridges_subset <- plot_multiset_gene_ridges(df = MeanVars,
                                                         gene_sets = gene_sets_to_ridge_subset,
                                                         gene_column = 'gene_ID',
                                                         value_column = 'HS_logCPM_Mean',
                                                         x_label = 'HS mean logCPM',
                                                         collapse_label = T)
multiset_gene_ridges_subset
ggsave(plot = multiset_gene_ridges_subset,filename = 'Subset_DEDV_MeanHS_Distributions.svg',width=2.5,height=2,dpi=300)

##############################################################
##### End) Save significant gene sets and workspace ##########
##############################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV")
write.csv(gamlss_NB_final, 'AllGenes_DE_DV.csv') 
write.csv(gene_set_numbers, 'DE_DV_GeneSet_Numbers.csv',row.names = F) 

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_GeneSets")

################## DEGs - different means #####################
write.csv(DEGsFinal, 'DEGs.csv') 
write.csv(DEGsIncreaseHS, 'DEGs_increased.csv')
write.csv(DEGsDecreaseHS, 'DEGs_decreased.csv')
# Large increases and decreases
write.csv(subset(DEGsIncreaseHS,logFC.cpm > 1),'DEGs_increased_large_FC.csv')
write.csv(subset(DEGsDecreaseHS,logFC.cpm < -1),'DEGs_decreased_large_FC.csv')

################## DVGs - different variance but not mean #####################
write.csv(DVGsFinal, 'DVGs.csv') 
write.csv(DVGsDecanalizedHS,'DVGs_decanalized.csv')
write.csv(DVGsCanalizedHS,'DVGs_canalized.csv')
write.csv(subset(DVGsDecanalizedHS,Difference.BCV > 0.5),'DVGs_decanalized_large_difference.csv')
write.csv(subset(DVGsCanalizedHS,Difference.BCV < -0.5),'DVGs_canalized_large_difference.csv')

################## NDE_NDV_Gs - same variabilities and means #####################
write.csv(NDE_NDV_Gs,'NDE_NDV_Gs.csv')

################## A(ll)DVGs - different variance regardless of mean change #####################
# All DEDVGs that either exhibit decanalization or canalization in HS
write.csv(ADVGsFinal, 'ADVGs.csv') 
write.csv(ADVGsDecanalizedHS,'ADVGs_decanalized.csv')
write.csv(ADVGsCanalizedHS,'ADVGs_canalized.csv')

################## DEDVGs - different variance and mean #####################
write.csv(DEDVGsFinal, 'DEDVGs.csv') 
write.csv(DEDVGsDecanalizedHS,'DEDVGs_decanalized.csv')
write.csv(DEDVGsCanalizedHS,'DEDVGs_canalized.csv')
write.csv(subset(DEDVGsDecanalizedHS,Difference.BCV > 0.5),'DEDVGs_decanalized_large_difference.csv')
# Note - I Flybased these 
# some of the most highly expressed genes of this set (10 CPM - 20 CPM) are pointed and frizzled
# These are the chief protein signals of the planar cell polarity pathway
# Implies a link to development of polar cells
# However, since there are many downstream players, it is difficult to draw big conclusions
# Imaging the polar cells of high-sugar vs non high-sugar flies might be doable...
write.csv(subset(DEDVGsCanalizedHS,Difference.BCV < 0.5),'DEDVGs_canalized_large_difference.csv')
write.csv(DEDVGsDecanalizedIncreasedHS,'DEDVGs_increased_decanalized.csv')
write.csv(DEDVGsDecanalizedDecreasedHS,'DEDVGs_decreased_decanalized.csv')
write.csv(DEDVGsCanalizedIncreasedHS,'DEDVGs_increased_canalized.csv')
write.csv(DEDVGsCanalizedDecreasedHS,'DEDVGs_decreased_canalized.csv')

# Save workspace
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV")
save.image(file='DE_DV_Main.RData') 


################## Post-script: weird genes breakdown #####################

# Plot logCPM (provided within the voom function) before and after batch effect correction
logCPM <- voom(raw.counts.DGElist[weird_genes,])
batch_effects <- model.matrix(~rnalib+rnaseq+egg+plate+well+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4)
treatment_effect <- model.matrix(~treatment)
logCPM.bc <- limma::removeBatchEffect(logCPM,covariates = batch_effects[,-1],design = treatment_effect) 
for (gene in weird_genes){
  par(mfcol=c(1,2))
  hist(logCPM$E[gene,Ctrl_flies],main = gene,xlab = "logCPM - Pre-Batch Correction")
  hist(logCPM$E[gene,HS_flies],add=T,col=rgb(0,1,1,0.5),main = gene)
  hist(logCPM.bc[gene,Ctrl_flies],main = gene,xlab = "logCPM - Post-Batch Correction")
  hist(logCPM.bc[gene,HS_flies],add=T,col=rgb(0,1,1,0.5),main = gene)}
# These are genes with a large range of RNA (-5 - 10). Some are bimodal before correction.
# Accounting for the batch effects (and not necessarily the SVs) is sufficient to turn the bimodal distributions normal

