# Transcriptome variability analyses of Adult Netherlands D.mel heads on normal vs high sugar diets
# Tests for differences in mean (DE) and variabilities (DV) between conditions

# Downsampling series
# Used to generate one % DE plot and one % DV plot
# Well is not used as a covariate here because often either only 1 or no samples of particular well are sampled, leading to non-convergence of the model

# Script by James Tan

# Last Updated: 10/9/25

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
library(gamlss)

# Obtain functions for mean-variance plots
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('Variability_Functions.R')

#####################
##### Datasets ######
#####################

# Metadata on all samples, e.g., conditions and batches
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",header = T)

# Surrogate variables, estimated in the investigating batch effect R scripts
# Based on manual PC inspection - 4 SVs for voom
tmm.voom.sv1.base <- read.table('TMM_Voom_sv1_10_wID.txt',row.names=1,header = T)
tmm.voom.sv2.base <- read.table('TMM_Voom_sv2_10_wID.txt',row.names=1,header = T)
tmm.voom.sv3.base <- read.table('TMM_Voom_sv3_10_wID.txt',row.names=1,header = T)
tmm.voom.sv4.base <- read.table('TMM_Voom_sv4_10_wID.txt',row.names=1,header = T)

# Expression - raw counts of filtered samples (see MakingGeneExpressionMatrix_head_HS&CTRL.R for filtering criteria)
raw.counts <- read.table("RawCounts_noY_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt",h=T,check.names = F)

# Eliminate bimodal genes 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
Bimodal_genes <- read.csv('Bimodal_Genes.csv',row.names = 1) 
Non_bimodal_genes <- setdiff(rownames(raw.counts),rownames(Bimodal_genes))
raw.counts <- raw.counts[Non_bimodal_genes,]

# Final dataframes to subset from 
cov <- tibble::column_to_rownames(head.info, "id")
RNASeq_flies <- cov[colnames(raw.counts),]
Ctrl_flies <- rownames(subset(RNASeq_flies,treatment==1))
HS_flies <- rownames(subset(RNASeq_flies,treatment==2))

################################
##### Prepare dataframes #######
################################

# Initiate the dataframes by including the full sample size numbers
downsampling_DE_df <- data.frame(n = 'Full', DEGs = 4382)
downsampling_DV_df <- data.frame(n = 'Full', DVGs = 7685)

##################################
##### Begin the subsampling ######
##################################

# What is the sample size?
downsample_sizes <- 500

# How many time do we randomly subsample? 
# We can code this as a set of starting seeds, one seed per subsample
downsampling_seeds <- 1:10

# Execute the double for loop to generate the necessary plot
for (sample_size in downsample_sizes){
  for (seed in downsampling_seeds){
    set.seed=seed
    Ctrl_flies_subsample <- sample(Ctrl_flies,size=sample_size,replace = F)
    HS_flies_subsample <- sample(HS_flies,size=sample_size,replace = F)
    RNASeq_flies_subsample <- c(Ctrl_flies_subsample,HS_flies_subsample)
    
    # Subset the counts and information to just the subsample
    raw.counts.sub <- raw.counts[,RNASeq_flies_subsample]
    cov2 <- cov[RNASeq_flies_subsample,]
    
    # Do the same for the surrogate variables
    tmm.voom.sv1 <- tmm.voom.sv1.base[RNASeq_flies_subsample,"x"]
    tmm.voom.sv2 <- tmm.voom.sv2.base[RNASeq_flies_subsample,"x"]
    tmm.voom.sv3 <- tmm.voom.sv3.base[RNASeq_flies_subsample,"x"]
    tmm.voom.sv4 <- tmm.voom.sv4.base[RNASeq_flies_subsample,"x"]
    
    # Merge expression and metadata 
    raw.counts.sub.t <- data.frame(t(raw.counts.sub)) 
    raw.counts.sub.t <- tibble::rownames_to_column(raw.counts.sub.t, "id")
    head.data <- merge(raw.counts.sub.t,head.info,by = 'id')
    head.data <- tibble::column_to_rownames(head.data, "id")
    
    # Preparing independent variable and batch variables
    cov2$treatment <- as.factor(cov2$treatment)
    cov2$RNAlibBatch <- as.factor(cov2$RNAlibBatch)
    cov2$RNAseqBatch <- as.factor(cov2$RNAseqBatch)
    cov2$egglayBatch<- as.factor(cov2$egglayBatch)
    cov2$platingBatch <- as.factor(cov2$platingBatch)
    #cov2$well <- as.factor(cov2$well)
    
    treatment <- cov2$treatment
    rnalib <- cov2$RNAlibBatch
    rnaseq <- cov2$RNAseqBatch
    egg <-cov2$egglayBatch
    plate <-cov2$platingBatch
    #well <-cov2$well
    
    # Sums to zero contrast for the batch effects
    contrasts(rnalib) <- contr.sum(levels(rnalib)) 
    contrasts(rnaseq) <- contr.sum(levels(rnaseq))
    contrasts(egg) <- contr.sum(levels(egg))
    contrasts(plate) <- contr.sum(levels(plate))
    #contrasts(well) <- contr.sum(levels(well))
    
    ##################################
    ##### Run GAMLSS DE and DV #######
    ##################################
    
    raw.counts.sub.DGElist <- DGEList(counts=raw.counts.sub)
    
    # Calculate TMM library size normalisation factors 
    raw.counts.sub.DGElist <- calcNormFactors(raw.counts.sub.DGElist, method="TMM")
    
    # Extract TMM-normalised library sizes to explicitly account for this in the model
    raw.counts.sub.DGElist$samples$offset <- log(raw.counts.sub.DGElist$samples$lib.size * raw.counts.sub.DGElist$samples$norm.factors)
    offset_values <- raw.counts.sub.DGElist$samples$offset
    
    # Full design matrix for model
    design <- model.matrix(~rnalib+rnaseq+egg+plate+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4+treatment)
    # Note: intercept represents the mean of control flies
    
    # GAMLSS Modelling
    
    # GAMLSS function
    gamlss_function <- function(gene,sigma.start.estimate,max_cycles) {
      y=raw.counts.sub.DGElist$counts[gene,]
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
        res$logFC.cpm = log(res$cpm.hs/res$cpm.ctrl)
        
        # GAMLSS log-likelihood ratio (LR) test for a significance of a treatnent effect on gene’s mean (CPM) counts.
        res$pval.cpm = pchisq(2*(logLik(m0)-logLik(m2)), df=m0$df.fit-m2$df.fit, lower=F)[[1]]
        
        # Write GAMLSS estimations of gene’s non-Poisson noise from the m0 model: bcv(μ)=√α.
        res$BCV.ctrl = sqrt(exp(m0$sigma.coefficients[[1]]))
        res$BCV.hs = sqrt(exp(m0$sigma.coefficients[[1]]+m0$sigma.coefficients[["treatment2"]]))
        res$sigma.ctrl = m0$sigma.coefficients[[1]]
        res$sigma.hs = m0$sigma.coefficients[[1]]+m0$sigma.coefficients[["treatment2"]]
        
        # Calculate log ratio for changes in cv(μ) between Ctrl and HS flies
        res$logFC.BCV = log(res$BCV.hs/res$BCV.ctrl)
        
        # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on non-Poisson noise.
        # padj.BCV – p value of LR test statistic: D_α=-2log⁡[L(μ_j,α_0  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m1 models.
        res$pval.BCV = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
      }
      res
    }
    
    # The default settings are 0.1 as a starting alpha
    # As well as 100 cycles maximum to find the best-fitting alpha
    gene_id <- 1:nrow(raw.counts.sub.DGElist)
    gamlss_NB <- lapply(gene_id, gamlss_function, sigma.start.estimate = 0.1, max_cycles=100)
    print(paste0("Seed ",seed," Complete"))
    
    # Merge the files
    gamlss_NB_2 <- do.call(rbind, gamlss_NB)
    rownames(gamlss_NB_2) <- rownames(raw.counts.sub.DGElist$counts)[gene_id]
    
    # Just keep all rows - do not account for weird genes or NA genes...
    gamlss_NB_final <- gamlss_NB_2[complete.cases(gamlss_NB_2),]
    nrow(gamlss_NB_final) 
    
    # Apply the BH correction
    gamlss_NB_final$padj.cpm <- p.adjust(gamlss_NB_final$pval.cpm, method = 'BH')
    gamlss_NB_final$padj.BCV <- p.adjust(gamlss_NB_final$pval.BCV, method = 'BH')
    
    ##############################################
    ##### 5) How many genes are DE and DV ########
    ##############################################
    
    # Save DE genes 
    DEGsFinal <- subset(gamlss_NB_final, padj.cpm < 0.05)
    DVGsFinal <- subset(gamlss_NB_final, padj.BCV < 0.05)
    
    # Add the number of DE genes and DV genes to the dataframe
    downsampling_DE_df<- rbind(data.frame(n = as.character(sample_size), DEGs = nrow(DEGsFinal)),downsampling_DE_df)
    downsampling_DV_df<- rbind(data.frame(n = as.character(sample_size), DVGs = nrow(DVGsFinal)),downsampling_DV_df)
  }
}

############################
##### Plot the results #####
############################

# Look at the filled dataframes
View(downsampling_DE_df) 
View(downsampling_DV_df) 

# Plot number of DEGs
DEG_number_by_subsample_size <- ggplot(data = downsampling_DE_df, aes(x=as.factor(n),y=DEGs)) + 
  geom_boxplot(fill='#4F8E4D',col='darkgrey')+
  geom_point(col='darkgrey')+
  theme_classic()+  
  ylab('Number of DEGs')+
  xlab('Sample size per condition')+  
  ylim(0,8763)+
  geom_hline(yintercept=8763, col="red")
DEG_number_by_subsample_size

DVG_number_by_subsample_size <- ggplot(data = downsampling_DV_df, aes(x=as.factor(n),y=DVGs)) + 
  geom_boxplot(fill='#611BB8',col='darkgrey')+
  geom_point(col='darkgrey')+
  theme_classic()+
  ylab('Number of DVGs')+
  xlab('Sample size per condition')+
  ylim(0,8763)+
  geom_hline(yintercept=8763, col="red")
DVG_number_by_subsample_size

# As a fraction of the full sample size
Fraction_DEG_by_subsample_size <- ggplot(data = downsampling_DE_df, aes(x=as.factor(n),y=100*DEGs/8763)) + 
  geom_boxplot(fill='#4F8E4D',col='darkgrey')+
  geom_point(col='darkgrey')+
  theme_classic()+  
  ylab('% of transcriptome DE')+
  xlab('Sample size per condition')+  
  ylim(0,100)
Fraction_DEG_by_subsample_size

Fraction_DVG_by_subsample_size <- ggplot(data = downsampling_DV_df, aes(x=as.factor(n),y=100*DVGs/8763)) + 
  geom_boxplot(fill='#611BB8',col='darkgrey')+
  geom_point(col='darkgrey')+
  theme_classic()+
  ylab('% of transcriptome DV')+
  xlab('Sample size per condition')+
  ylim(0,100)
Fraction_DVG_by_subsample_size

##################################
##### End) Save the workspace ####
##################################

# Save workspace
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_Downsampled/")
save.image(file='DE_DV_Downsampled_Series_Small_500.RData')
