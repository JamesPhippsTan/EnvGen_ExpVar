# Transcriptome variability analyses of Adult Netherlands D.mel heads on normal vs high sugar diets
# Tests the robustness of variability estimate to down-sampling

# Script by James Tan

# Last Updated: 12/9/25

################################
##### Packages at Setup ########
################################

rm(list = ls())

library(DESeq2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(car)
library(matrixStats)
library(ggplot2)
library(ggrepel)

# Working Directory 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD")

# Load saved script environment
load(file='Downsampling.RData')

# Functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code\\")
source('Variability_Functions.R')
source('Quantile_Functions.R')

# Function to downsample the expression matrix and compute variance
compute_downsampled_variabilities <- function(expression_matrix, sample_sizes, n_iterations) {
  downsampled_variabilities <- list()
  
  for (size in sample_sizes) {
    message("Processing sample size: ", size)
    
    # Matrix to store variance estimates for each iteration
    var_matrix <- matrix(nrow = nrow(expression_matrix), ncol = n_iterations)
    rownames(var_matrix) <- rownames(expression_matrix)
    
    for (i in 1:n_iterations) {
      set.seed(i)
      sampled_indices <- sample(ncol(expression_matrix), size = size, replace = F)
      var_matrix[, i] <- rowMads(expression_matrix[, sampled_indices],constant=1)
    }
    
    downsampled_variabilities[[paste0("Var_", size)]] <- var_matrix
  }
  
  return(downsampled_variabilities)
}

# Function to compute percent difference between full variance and subsampled variance
compute_percent_diff <- function(subsample_matrix, full_var) {
  100 * rowMeans(subsample_matrix - full_var) / full_var
}

# Function to initialize data frames dynamically
initialize_df <- function(mean_col, var_col, row_names, sample_sizes) {
  df <- data.frame(
    Mean = mean_col,
    Var = var_col,
    MeanQuantile = NA,
    VarQuantile = NA,
    row.names = row_names
  )
  
  for (size in sample_sizes) {
    df[[paste0("PercentDiff", size)]] <- NA
  }
  
  return(df)
}

# Function to compute percent differences for all sample sizes
compute_percent_differences <- function(var_list, full_var_df, sample_sizes) {
  for (size in sample_sizes) {
    var_name <- paste0("Var_", size)
    full_var_df[[paste0("PercentDiff", size)]] <- compute_percent_diff(
      var_list[[var_name]][rownames(full_var_df), ], full_var_df$Var
    )
  }
  return(full_var_df)
}


######################
##### Datasets #######
######################

# Metadata on all samples, e.g., conditions at batches
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)

# Surrogate variables, estimated in the SurrogateVariables R scripts
vst.sv1 <- read.table('VST_sv1_10.txt')
vst.sv1 <- vst.sv1$x
vst.sv2 <- read.table('VST_sv2_10.txt')
vst.sv2 <- vst.sv2$x
vst.sv3 <- read.table('VST_sv3_10.txt')
vst.sv3 <- vst.sv3$x

# Expression - raw counts of filtered samples (see MakingGeneExpressionMatrix_head_HS&CTRL.R for filtering criteria)
raw.counts <- read.table("RawCounts_noY_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt",h=T,check.names = F)

# Eliminate bimodal genes
Bimodal_genes <- read.csv('Bimodal_Genes.csv',row.names = 1) 
Non_bimodal_genes <- setdiff(rownames(raw.counts),rownames(Bimodal_genes))
raw.counts <- raw.counts[Non_bimodal_genes,]

# Only metadata on samples with counts
cov <- tibble::column_to_rownames(head.info, "id")
cov2 <- cov[c(colnames(raw.counts)),]

# Merge expression at metadata 
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

#####################
##### VST ###########
#####################

# Apply variance-stabilizing transformation
data <- DESeqDataSetFromMatrix(countData = raw.counts,
                               colData = cov2,
                               design =  ~treatment) 
vst.data = vst(data, blind=F) 

# Covariate matrices
VST.covariates <- model.matrix(~rnalib+rnaseq+egg+plate+well+vst.sv1+vst.sv2+vst.sv3)

# Correct counts for batch effects but not for treatment
vst.counts <- assay(vst.data)
vst.counts.bc  <- limma::removeBatchEffect(vst.counts, covariates = VST.covariates[,-1],  design = nullmodel) 

############################################################
##### 1) Estimating variability with all the samples #######
############################################################
# For this, we need to subset the raw expression data 

# Variance for each condition
conditions <- factor(t(head.data$treatment))
conditionsLevel <- levels(conditions)
Ctrl_flies <- c(which(conditions==conditionsLevel[1]))
HS_flies <- c(which(conditions==conditionsLevel[2]))

# Define sample sizes for downsampling
sample_sizes <- c(10, 50, 100, 200, 300, 500, 900)

# Define how many times to run downsampling
n_subsamples <- 100

# Split datasets
Condition_index_list = list(Ctrl=Ctrl_flies,HS=HS_flies)

for (condition_name in names(Condition_index_list)){
  
  condition <- Condition_index_list[[condition_name]]
    
  # Expression matrix (Replace with actual data)
  expression_matrix <- vst.counts.bc[,condition]
  
  # Compute full variance
  full_var <- rowMads(expression_matrix, constant = 1)
  full_mean <- rowMeans(expression_matrix)
  
  # Compute downsampled variabilities
  downsampled_variabilities <- compute_downsampled_variabilities(expression_matrix, sample_sizes,n_iterations = n_subsamples)
  
  # Initialize data frames
  VarDiffPercent <- initialize_df(full_mean, full_var, rownames(expression_matrix), sample_sizes)
  
  # Compute percent differences
  VarDiffPercent <- compute_percent_differences(downsampled_variabilities, VarDiffPercent, sample_sizes)
  
  # Compute quantiles for mean and variance
  MeanQuantileName <- paste0("Mean_Quantile_",condition_name)
  VarQuantileName <- paste0("Var_Quantile_",condition_name)
  VarDiffPercent <- make_quantiles(VarDiffPercent, VarDiffPercent$Mean, n_quantiles = 10, name_of_quantiles = MeanQuantileName)
  VarDiffPercent <- make_quantiles(VarDiffPercent, VarDiffPercent$Var, n_quantiles = 10, name_of_quantiles = VarQuantileName)
  
  # Convert to long format for plotting
  VarDiffPercent_p <- pivot_longer(
    data = VarDiffPercent,
    cols = starts_with("PercentDiff"),
    names_to = "Subsample_size",
    values_to = "Percent_Difference_Variability_Estimate"
  )
  
  # Convert Subsample_size to factor and rename
  VarDiffPercent_p$Subsample_size <- factor(VarDiffPercent_p$Subsample_size, levels = paste0("PercentDiff", sample_sizes), labels = as.character(sample_sizes))
  
  # Plot by mean quantile
  Subsample_VarDiff_byMean <- ggplot(VarDiffPercent_p, aes(y = Percent_Difference_Variability_Estimate, x = as.factor(!!sym(MeanQuantileName)), fill = Subsample_size)) +
    geom_boxplot() +
    theme_classic()+
    geom_hline(yintercept = c(5, -5), color = 'red') +
    ylab(paste0('Average % Difference from true variability - ',condition_name)) +
    xlab('Mean Quantile (876 genes)') +
    ggtitle("100 random subsamples with replacement")
  print(Subsample_VarDiff_byMean)
  ggsave(plot = Subsample_VarDiff_byMean,filename = paste0(condition_name,"_Subsample_VarDiff_byMean.png"),
         units = "in",
         width = 4.68, height = 4.68)

  # Plot by variability quantile
  Subsample_VarDiff_byVar <- ggplot(VarDiffPercent_p, aes(y = Percent_Difference_Variability_Estimate, x = as.factor(!!sym(VarQuantileName)), fill = Subsample_size)) +
    geom_boxplot() +
    theme_classic()+
    geom_hline(yintercept = c(5, -5), color = 'red') +
    ylab(paste0('Average % difference from true variability - ',condition_name)) +
    xlab('Variability Quantile (876 genes)') 
  print(Subsample_VarDiff_byVar)
  ggsave(plot = Subsample_VarDiff_byVar,filename = paste0(condition_name,"_Subsample_VarDiff_byVar.png"),
         units = "in",
         width = 4, height = 5)
  

  # Overall boxplot by sample size
  Subsample_VarDiff_bySampleSize <- ggplot(VarDiffPercent_p, aes(y = Percent_Difference_Variability_Estimate, x = Subsample_size)) +
    geom_boxplot() +
    theme_classic()+
    ylab(paste0('Average % difference from true variability (',condition_name,')')) +
    xlab('Subsample Size') 
  ggsave(dpi = 300,plot = Subsample_VarDiff_bySampleSize,
         filename = paste0(condition_name,"_Subsample_VarDiff_bySampleSize.svg"),
         units = "in",
         width = 4, height = 5)
}

# Save workspace
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD")
save.image(file='Downsampling.RData') 

