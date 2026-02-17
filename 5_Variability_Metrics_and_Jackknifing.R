# Transcriptome variance analyses of Adult Netherlands D.mel heads on normal vs high sugar diets
# Variance metrics and effect of adjustment
# Used to infer directionality of variance changes 
# Also used to determine variance rankings to be correlated with other genic features

# The best variance metric should
# 1) Have minimal correlation with mean - tested by correlation tests
# 2) Be minimally impacted by outliers - tested by the range of percent and rank differences between jackknifed metrics and full dataset metrics

# Script by James Tan
# Last Updated: 10/2/25

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
library(DescTools)
library(bootstrap)


# Obtain functions used to calculate adjusted SD and bootstraps
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('Variability_Functions.R')

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
#load(file='Variability_Metrics_and_Jackknifing.RData')

#####################
##### Datasets ######
#####################

# Metadata on all samples, e.g., conditions and batches
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)

# Surrogate variables, estimated in the NvsHS_SurrogateVariables R scripts
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

contrasts(rnalib) <- contr.sum(levels(rnalib)) # sums to zero contrast
contrasts(rnaseq) <- contr.sum(levels(rnaseq))
contrasts(egg) <- contr.sum(levels(egg))
contrasts(plate) <- contr.sum(levels(plate))
contrasts(well) <- contr.sum(levels(well))


# Defining indices to subset the total dataset into Ctrl and HS downstream
conditions <- factor(t(head.data$treatment))
conditionsLevel <- levels(conditions)
Ctrl_flies <- c(which(conditions==conditionsLevel[1]))
HS_flies <- c(which(conditions==conditionsLevel[2]))



# nullmodel - the biological factor of interest
nullmodel <- model.matrix(~treatment)

######################
######### VST ########
######################

data <- DESeqDataSetFromMatrix(countData = raw.counts,
                               colData = cov2,
                               design = ~treatment) 

# Apply VST()
vst.data = vst(data, blind=F, fitType = 'parametric') 
vst.counts <- assay(vst.data)

# Remove batch effects 
VST.covariates <- model.matrix(~rnalib+rnaseq+egg+plate+well+vst.sv1+vst.sv2+vst.sv3)
vst.counts.bc  <- limma::removeBatchEffect(vst.counts, covariates = VST.covariates[,-1],  design = nullmodel) 

##################################################
######### What variance metric to use? ###########
##################################################

# Variance vs mean
Variability_mean_plots(exp_matrix = vst.counts.bc,samples = Ctrl_flies,var_func = rowVars,mean_func = rowMeans,correlation_method = "spearman",log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")
Variability_mean_plots(exp_matrix = vst.counts.bc,samples = HS_flies,var_func = rowVars,mean_func = rowMeans,correlation_method = "spearman",log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")

# SD vs mean
Variability_mean_plots(exp_matrix = vst.counts.bc,samples = Ctrl_flies,var_func = rowSds,mean_func = rowMeans,correlation_method = "spearman",log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "SD")
Variability_mean_plots(exp_matrix = vst.counts.bc,samples = HS_flies,var_func = rowSds,mean_func = rowMeans,correlation_method = "spearman",log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "SD")

# IQRs vs mean
Variability_mean_plots(exp_matrix = vst.counts.bc,samples = Ctrl_flies,var_func = rowIQRs,mean_func = rowMeans,correlation_method = "spearman",log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "IQR")
Variability_mean_plots(exp_matrix = vst.counts.bc,samples = HS_flies,var_func = rowIQRs,mean_func = rowMeans,correlation_method = "spearman",log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "IQR")

# MADs vs median
Variability_mean_plots(exp_matrix = vst.counts.bc,samples = Ctrl_flies,var_func = rowMads,mean_func = rowMedians,correlation_method = "spearman",log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Median absolute deviation")
Variability_mean_plots(exp_matrix = vst.counts.bc,samples = HS_flies,var_func = rowMads,mean_func = rowMedians,correlation_method = "spearman",log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Median absolute deviation")

# The MAD has the lowest spearman correlations

###############################################################
######### Jackknifing effect on the metrics  ##################
###############################################################

jackknifed_Means_Ctrl <- jackknife_tables(vst.counts.bc,samples = Ctrl_flies,func = rowMeans)
jackknifed_Variances_Ctrl <- jackknife_tables(vst.counts.bc,samples = Ctrl_flies,func = rowVars)
jackknifed_SD_Ctrl <- jackknife_tables(vst.counts.bc,samples = Ctrl_flies,func = rowSds)
jackknifed_IQR_Ctrl <- jackknife_tables(vst.counts.bc,samples = Ctrl_flies,func = rowIQRs)
jackknifed_MAD_Ctrl <- jackknife_tables(vst.counts.bc,samples = Ctrl_flies,func = rowMads)

jackknifed_Means_HS <- jackknife_tables(vst.counts.bc,samples = HS_flies,func = rowMeans)
jackknifed_Variances_HS <- jackknife_tables(vst.counts.bc,samples = HS_flies,func = rowVars)
jackknifed_SD_HS <- jackknife_tables(vst.counts.bc,samples = HS_flies,func = rowSds)
jackknifed_IQR_HS <- jackknife_tables(vst.counts.bc,samples = HS_flies,func = rowIQRs)
jackknifed_MAD_HS <- jackknife_tables(vst.counts.bc,samples = HS_flies,func = rowMads)

# For each jackknifed object, extract the rank difference ranges 
jackknifed_stat_list <- mget(ls(pattern = "^jackknifed_"))

# Rank differences 
Rank_differences_range <- do.call(rbind, lapply(names(jackknifed_stat_list), function(name) {
  data.frame(name = gsub(x = name, pattern = "jackknifed_",replacement = ""), 
             max_rank_decrease = jackknifed_stat_list[[name]]$rank_differences_range[1],
             max_rank_increase = jackknifed_stat_list[[name]]$rank_differences_range[2])
}))

# Percent differences
Percent_differences_range <- do.call(rbind, lapply(names(jackknifed_stat_list), function(name) {
  data.frame(name = gsub(x = name, pattern = "jackknifed_",replacement = ""), 
             max_percent_decrease = jackknifed_stat_list[[name]]$percent_differences_range[1],
             max_percent_increase = jackknifed_stat_list[[name]]$percent_differences_range[2])
}))

# There is a smaller range for the MAD - therefore I will use this


####################################
######### End) Save files  #########
####################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")

# Save the range of jackknife estimates
write.csv(Rank_differences_range,"Jackknife_rank_differences.csv")
write.csv(Percent_differences_range,"Jackknife_percent_differences.csv")

# Save workspace
save.image(file='Variability_Metrics_and_Jackknifing.RData') 

