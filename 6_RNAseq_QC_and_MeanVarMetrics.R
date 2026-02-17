# Check various quality metrics of the final RNAseq dataset
# Try different transformations and metrics of mean and variability to identify the most suitable ones
# Afterwards, create a dataframes containing the means and variabilities
# Script by James Tan

# Last Updated: 29/7/25

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
library(gridExtra)
library(patchwork)

# Load file
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
load(file='Decide_MeanVar_Metrics.RData')

# Functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('Variability_Functions.R')

#######################
##### Datasets ########
#######################

# Metadata on all samples, e.g., conditions and batches
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)

# Surrogate variables, estimated in the SurrogateVariables R scripts
vst.sv1 <- read.table('VST_sv1_10.txt')
vst.sv1 <- vst.sv1$x
vst.sv2 <- read.table('VST_sv2_10.txt')
vst.sv2 <- vst.sv2$x
vst.sv3 <- read.table('VST_sv3_10.txt')
vst.sv3 <- vst.sv3$x

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
Bimodal_genes <- read.csv('Bimodal_Genes.csv',row.names = 1) 
Non_bimodal_genes <- setdiff(rownames(raw.counts),rownames(Bimodal_genes))
raw.counts <- raw.counts[Non_bimodal_genes,]

###############################
##### Process datasets ########
###############################

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

# Preparing treatment and variables
cov2$treatment <- as.factor(cov2$treatment)
cov2$well <- as.factor(cov2$well)
cov2$plate <- as.factor(cov2$plate)
cov2$RNAlibBatch <- as.factor(cov2$RNAlibBatch)
cov2$RNAseqBatch <- as.factor(cov2$RNAseqBatch)
cov2$egglayBatch<- as.factor(cov2$egglayBatch)
cov2$platingBatch <- as.factor(cov2$platingBatch)
cov2$eclosionBatch <- as.factor(cov2$eclosionBatch)
treatment <- cov2$treatment
well <- cov2$well
rnalib <- cov2$RNAlibBatch
rnaseq <- cov2$RNAseqBatch
egg <-cov2$egglayBatch
plate <- cov2$platingBatch
contrasts(well) <- contr.sum(levels(well))
contrasts(rnalib) <- contr.sum(levels(rnalib))
contrasts(rnaseq) <- contr.sum(levels(rnaseq))
contrasts(egg) <- contr.sum(levels(egg))
contrasts(plate) <- contr.sum(levels(plate))
batch_effects <- model.matrix(~rnalib+rnaseq+egg+plate+well+vst.sv1+vst.sv2+vst.sv3)
treatment_model <- model.matrix(~treatment)

########################
##### RNAseq QC ########
########################

# Set up 
Condition_df <- c(rep("Control", length(Ctrl_flies)), 
                  rep("HS", length(HS_flies)))

# Prefiltered read number per fly of expressed genes
prefiltered.raw.counts.Ctrl <- head.data[Ctrl_flies,'total_reads']/1000000
prefiltered.raw.counts.HS <- head.data[HS_flies,'total_reads']/1000000
prefiltered.raw.counts.df <- data.frame(Counts_per_fly = c(prefiltered.raw.counts.Ctrl,prefiltered.raw.counts.HS),
                            Condition = Condition_df)
prefiltered.raw.counts.plot <- ggplot(prefiltered.raw.counts.df, aes(y=Counts_per_fly,x = Condition,fill=Condition)) + 
  geom_boxplot(col='darkgrey') +
  ylab(expression('Total read counts per fly / x'~10^6))+
  xlab('Condition')+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  theme_classic()  +
  scale_fill_manual(values = c("Control" = "#C6B49E", "HS" = "#DF9F65"))
prefiltered.raw.counts.plot

# Raw read number per fly of expressed genes - after filtering
raw.counts.Ctrl <- colSums(raw.counts[,Ctrl_flies])/1000000
raw.counts.HS <- colSums(raw.counts[,HS_flies])/1000000
raw.counts.df <- data.frame(Counts_per_fly = c(raw.counts.Ctrl,raw.counts.HS),
                            Condition = Condition_df)
raw.counts.plot <- ggplot(raw.counts.df, aes(y=Counts_per_fly,x = Condition,fill=Condition)) + 
  geom_boxplot(col='darkgrey') +
  ylab(expression('Filtered read counts per fly / x'~10^6))+
  xlab('Condition')+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  theme_classic()  +
  scale_fill_manual(values = c("Control" = "#C6B49E", "HS" = "#DF9F65"))
raw.counts.plot


# Mean normalised expression per fly
raw.counts.list <- DGEList(counts=raw.counts)
TMM.counts <- calcNormFactors(raw.counts.list, method = 'TMM') 
voom.object <- voom(TMM.counts)
TMM.Voom.covariates <- model.matrix(~rnalib+rnaseq+egg+plate+well+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4)
voom.counts.bc <- limma::removeBatchEffect(voom.object, 
                                           covariates = TMM.Voom.covariates[,-1],  
                                           design = treatment_model) 
voom.counts.bc <- as.data.frame(voom.counts.bc)
voom.counts.Ctrl.fly <- colMeans(voom.counts.bc[,Ctrl_flies])
voom.counts.HS.fly <- colMeans(voom.counts.bc[,HS_flies])
voom.counts.df.fly <- data.frame(Counts_per_fly = c(voom.counts.Ctrl.fly,voom.counts.HS.fly),
                            Condition = Condition_df)
voom.counts.plot.fly <- ggplot(voom.counts.df.fly, aes(y=Counts_per_fly,x = Condition,fill=Condition)) + 
  geom_boxplot(col='darkgrey') +
  ylab(expression('Average voom expression per fly'))+
  xlab('Condition')+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  theme_classic()  +
  scale_fill_manual(values = c("Control" = "#C6B49E", "HS" = "#DF9F65"))
voom.counts.plot.fly

# Mean normalised expression per gene
voom.counts.Ctrl <- rowMeans(voom.counts.bc[,Ctrl_flies])
voom.counts.HS <- rowMeans(voom.counts.bc[,HS_flies])
Condition_df_gene <- c(rep("Control", nrow(voom.counts.bc)), 
                  rep("HS", nrow(voom.counts.bc)))
voom.counts.df <- data.frame(Counts_per_fly = c(voom.counts.Ctrl,voom.counts.HS),
                             Condition = Condition_df_gene)
voom.counts.plot <- ggplot(voom.counts.df, aes(y=Counts_per_fly,x = Condition,fill=Condition)) + 
  geom_boxplot(col='darkgrey') +
  ylab(expression('Voom expression per gene'))+
  xlab('Condition')+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  theme_classic()  +
  scale_fill_manual(values = c("Control" = "#C6B49E", "HS" = "#DF9F65"))
voom.counts.plot

QC_plots <-  prefiltered.raw.counts.plot+ theme(legend.position = "none") + 
  raw.counts.plot+ theme(legend.position = "none") + 
  voom.counts.plot.fly + theme(legend.position = "none")+ 
  voom.counts.plot+ theme(legend.position = "none")+ 
  plot_layout(ncol = 2)
QC_plots

###############################################
##### Different variabilites ##################
###############################################

# Now that we know the data is useable, let us try different metrics of means and variabilities

# Mean-variance plots before any transformation
Ctrl_MeanVar_raw_plots <- Variability_mean_plots(exp_matrix = raw.counts,samples = Ctrl_flies,var_func = rowMads,mean_func = rowMeans,log_mean = T,log_var = T,mean_name = "Mean log10(counts) (Ctrl)",var_name = "log10(Variability) (MAD)")
Ctrl_MeanVar_raw_plots_full <- Ctrl_MeanVar_raw_plots[[1]]
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = Ctrl_MeanVar_raw_plots_full,filename = 'Ctrl_MeanVar_raw_plots_full.svg',dpi=300,width=2.24,height = 2.24)

HS_MeanVar_raw_plots <- Variability_mean_plots(exp_matrix = raw.counts,samples = HS_flies,var_func = rowMads,mean_func = rowMeans,log_mean = T,log_var = T,mean_name = "Mean log10(counts) (HS)",var_name = "log10(Variability) (MAD)")
HS_MeanVar_raw_plots_full <- HS_MeanVar_raw_plots[[1]]
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = HS_MeanVar_raw_plots_full,filename = 'HS_MeanVar_raw_plots_full.svg',dpi=300,width=2.24,height = 2.24)


###################################################
##### (a) Different VST transformations ###########
###################################################

# Parametric (default) fit
data <- DESeqDataSetFromMatrix(countData = raw.counts,
                               colData = cov2,
                               design =  ~treatment) 
vst.data = vst(data, blind=F,fitType = "parametric")
vst.countsm <- as.matrix(assay(vst.data))
counts.adj_parametric <- limma::removeBatchEffect(vst.countsm, covariates = batch_effects[,-1],  design = treatment_model) 
Variability_mean_plots(exp_matrix = counts.adj_parametric,samples = Ctrl_flies,var_func = rowVars,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")
Variability_mean_plots(exp_matrix = counts.adj_parametric,samples = HS_flies,var_func = rowVars,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")
Variability_mean_plots(exp_matrix = counts.adj_parametric,samples = Ctrl_flies,var_func = rowMads,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Median logCPM",var_name = "MAD")
Variability_mean_plots(exp_matrix = counts.adj_parametric,samples = HS_flies,var_func = rowMads,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Median logCPM",var_name = "MAD")

# Local fit
vst.data = vst(data, blind=F,fitType = "local") 
vst.countsm <- as.matrix(assay(vst.data))
counts.adj_local <- limma::removeBatchEffect(vst.countsm, covariates = batch_effects[,-1],  design = treatment_model) 
Variability_mean_plots(exp_matrix = counts.adj_local,samples = Ctrl_flies,var_func = rowVars,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")
Variability_mean_plots(exp_matrix = counts.adj_local,samples = HS_flies,var_func = rowVars,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")
# Similar to base parametric VST

# Mean fit
vst.data = vst(data, blind=F,fitType = "mean") 
vst.countsm <- as.matrix(assay(vst.data))
counts.adj_meanfit <- limma::removeBatchEffect(vst.countsm, covariates = batch_effects[,-1],  design = treatment_model) 
Variability_mean_plots(exp_matrix = counts.adj_meanfit,samples = Ctrl_flies,var_func = rowVars,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")
Variability_mean_plots(exp_matrix = counts.adj_meanfit,samples = HS_flies,var_func = rowVars,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")
# Virtually the same plot as the BCV against the logCPM

# Parametric-fitting Ctrl and HS separately
# Estimation with knowledge of the design (blind=F), as above
objectblind_F <- DESeqDataSetFromMatrix(countData = raw.counts,
                                   colData = cov2,
                                   design =  ~treatment) 
objectblind_F <- estimateSizeFactors(objectblind_F)
objectblind_F <- estimateDispersionsGeneEst(objectblind_F)
design(objectblind_F)
objectblind_F_fit <- estimateDispersionsFit(objectblind_F, quiet=TRUE, fitType = "parametric")
plotDispEsts(objectblind_F_fit)
environment(objectblind_F_fit@dispersionFunction)[["coefs"]]
#asymptDisp  extraPois 
#0.1524928 23.9270910 

# Estimation with just the Ctrl dataset
object_Ctrl <- DESeqDataSetFromMatrix(countData = raw.counts[,Ctrl_flies],
                                      colData = cov2[Ctrl_flies,],design = ~1) 
object_Ctrl <- estimateSizeFactors(object_Ctrl)
object_Ctrl <- estimateDispersionsGeneEst(object_Ctrl)
object_Ctrl_fit <- estimateDispersionsFit(object_Ctrl, quiet=TRUE, fitType = "parametric")
plotDispEsts(object_Ctrl_fit)
environment(object_Ctrl_fit@dispersionFunction)[["coefs"]]
# asymptDisp  extraPois 
# 0.1193211 12.9929837 
# A lot lower than the average and the HS
# Estimation with just the HS dataset
object_HS <- DESeqDataSetFromMatrix(countData = raw.counts[,HS_flies],
                                    colData = cov2[HS_flies,],design = ~1) 
object_HS <- estimateSizeFactors(object_HS)
object_HS <- estimateDispersionsGeneEst(object_HS)
object_HS_fit <- estimateDispersionsFit(object_HS, quiet=TRUE, fitType = "parametric")
plotDispEsts(object_HS_fit)
environment(object_HS_fit@dispersionFunction)[["coefs"]]
# asymptDisp  extraPois 
# 0.1989378 28.6237915 
# A lot higher than the average and the Ctrl
# Transform separately and then re-collate the dataframe
Ctrl_vsd <- getVarianceStabilizedData(object_Ctrl_fit)
HS_vsd <- getVarianceStabilizedData(object_HS_fit)
Ctrl_HS_indiv_VST <- cbind(Ctrl_vsd,HS_vsd)
Ctrl_HS_indiv_VST <- Ctrl_HS_indiv_VST[,colnames(vst.countsm)]
counts.adj_indiv_VST <- limma::removeBatchEffect(Ctrl_HS_indiv_VST, covariates = batch_effects[,-1],  design = treatment_model) 
Variability_mean_plots(exp_matrix = counts.adj_indiv_VST,samples = Ctrl_flies,var_func = rowVars,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")
Variability_mean_plots(exp_matrix = counts.adj_indiv_VST,samples = HS_flies,var_func = rowVars,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM",var_name = "Variance")
# The trend is still present - even using the condition-specific stabilization does not eliminate this

# Final decided metric for plots - VST 
Ctrl_MeanVar_VST_plots <- Variability_mean_plots(exp_matrix = counts.adj_parametric,samples = Ctrl_flies,var_func = rowMads,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean of VST counts (Ctrl)",var_name = "log10(Variability) (MAD)")
Ctrl_MeanVar_VST_plots_full <- Ctrl_MeanVar_VST_plots[[1]]
Ctrl_MeanVar_VST_plots_full
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = Ctrl_MeanVar_VST_plots_full,filename = 'Ctrl_MeanVar_VST_plots_full.svg',dpi=300,width=2.24,height = 2.24)

HS_MeanVar_VST_plots <- Variability_mean_plots(exp_matrix = counts.adj_parametric,samples = HS_flies,var_func = rowMads,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean of VST counts (HS)",var_name = "log10(Variability) (MAD)")
HS_MeanVar_VST_plots_full <- HS_MeanVar_VST_plots[[1]]
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = HS_MeanVar_VST_plots_full,filename = 'HS_MeanVar_VST_plots_full.svg',dpi=300,width=2.24,height = 2.24)

########################################################
##### (b) Adjusted variabilities based on logCPM #######
########################################################

# Final decided metric for plots - logCPM 
Ctrl_MeanVar_logCPM_plots <- Variability_mean_plots(exp_matrix = voom.counts.bc,samples = Ctrl_flies,var_func = rowMads,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM (Ctrl)",var_name = "log10(Variability) (MAD)")
Ctrl_MeanVar_logCPM_plots_full <- Ctrl_MeanVar_logCPM_plots[[1]]
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = Ctrl_MeanVar_logCPM_plots_full,filename = 'Ctrl_MeanVar_logCPM_plots_full.svg',dpi=300,width=2.24,height = 2.24)

HS_MeanVar_logCPM_plots <- Variability_mean_plots(exp_matrix = voom.counts.bc,samples = HS_flies,var_func = rowMads,mean_func = rowMeans,log_mean = F,log_var = T,mean_name = "Mean logCPM (HS)",var_name = "log10(Variability) (MAD)")
HS_MeanVar_logCPM_plots_full <- HS_MeanVar_logCPM_plots[[1]]
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = HS_MeanVar_logCPM_plots_full,filename = 'HS_MeanVar_logCPM_plots_full.svg',dpi=300,width=2.24,height = 2.24)

# Prepare table to store mean-variance and between-metric correlations
VST_Metric_Correlation_Table <- data.frame(
  Condition = character(),
  VST_MAD_MeanVar_Spearman_Coeff = numeric(),
  VST_MAD_MeanVar_Spearman_pval = numeric(),
  Adj_Var_Polynomial = numeric(),
  Adj_Var_MeanVar_Spearman_Coeff = numeric(),
  Adj_Var_MeanVar_Spearman_pval = numeric(),
  LOESS_Var_MeanVar_Spearman_Coeff = numeric(),
  LOESS_Var_MeanVar_Spearman_pval = numeric(),
  VST_MAD_Adj_Var_Spearman_Coeff = numeric(),
  VST_MAD_Adj_Var_Spearman_pval = numeric(),
  VST_MAD_LOESS_Var_Spearman_Coeff = numeric(),
  VST_MAD_LOESS_Var_Spearman_pval = numeric(),
  stringsAsFactors = FALSE
)

# Prepare table to store variance and mean metrics
Mean_Var_Table <- data.frame(
  gene_ID = rownames(voom.counts.bc)
)

# Get metrics
Condition_index_list = list(Ctrl=Ctrl_flies,HS=HS_flies)
for (condition_name in names(Condition_index_list)){
  
  condition <- Condition_index_list[[condition_name]]
  
  # VST() MADs
  VST_MADs <- rowMads(counts.adj_parametric[,condition],constant=1)
  VST_Means <- rowMeans(counts.adj_parametric[,condition])
  VST_MAD_MeanVar <- cor.test(VST_MADs,VST_Means,method='spearman')
  
  # Global and condition-specific variabilities and means
  Global_Vars <- rowVars(voom.counts.bc)
  Global_Means <- rowMeans(voom.counts.bc)
  Vars <- rowVars(voom.counts.bc[,condition])
  Means <- rowMeans(voom.counts.bc[,condition])
  
  # Adjusted variance - Polynomial regression method
  adjusted_var <- polynomial_adjust_var(global_means = Global_Means,
                                        global_vars = Global_Vars,
                                        specific_vars = Vars)
  plot(Means,log(adjusted_var$Vars))
  plot(rank(Means),rank(adjusted_var$Vars))
  Adj_Var_MeanVar <- cor.test(Means,adjusted_var$Vars,method='spearman')

  # LOESS Method
  loess_var = loess(log(Vars) ~ Means, degree = 1, span = 0.7)
  loess_vars = resid(loess_var)
  plot(Means,loess_vars)
  plot(rank(Means),rank(loess_vars))
  LOESS_Var_MeanVar <- cor.test(Means,loess_vars,method='spearman')
  
  # Plot against other types of VSTs
  plot(rank(VST_MADs),rank(adjusted_var$Vars))
  plot(rank(VST_MADs),rank(loess_vars))
  
  # Get the correlations
  VST_MAD_LOESS_Var <- cor.test(VST_MADs,loess_vars,method='spearman')
  VST_MAD_Adj_Var <- cor.test(VST_MADs,adjusted_var$Vars,method='spearman')
  
  # Save the correlations between the different metrics
  VST_Metric_Correlation_Table <- rbind(VST_Metric_Correlation_Table,data.frame(
    Condition = condition_name,
    VST_MAD_MeanVar_Spearman_Coeff = format_statistic(VST_MAD_MeanVar$estimate),
    VST_MAD_MeanVar_Spearman_pval = format_statistic(VST_MAD_MeanVar$p.value),
    Adj_Var_Polynomial = format_statistic(adjusted_var$polynomial),
    Adj_Var_MeanVar_Spearman_Coeff = format_statistic(Adj_Var_MeanVar$estimate),
    Adj_Var_MeanVar_Spearman_pval = format_statistic(Adj_Var_MeanVar$p.value),
    LOESS_Var_MeanVar_Spearman_Coeff = format_statistic(LOESS_Var_MeanVar$estimate),
    LOESS_Var_MeanVar_Spearman_pval = format_statistic(LOESS_Var_MeanVar$p.value),
    VST_MAD_Adj_Var_Spearman_Coeff = format_statistic(VST_MAD_Adj_Var$estimate),
    VST_MAD_Adj_Var_Spearman_pval = format_statistic(VST_MAD_Adj_Var$p.value),
    VST_MAD_LOESS_Var_Spearman_Coeff = format_statistic(VST_MAD_LOESS_Var$estimate),
    VST_MAD_LOESS_Var_Spearman_pval = format_statistic(VST_MAD_LOESS_Var$p.value),
    stringsAsFactors = FALSE
  ))
  
  # Add the metrics to the table
  Metrics = data.frame(
      logCPM_Variance = Vars,
      logCPM_Mean = Means,
      VST_MAD = VST_MADs,
      VST_Mean = VST_Means,
      Adj_Var = adjusted_var$Vars,
      LOESS_Var = loess_vars
  )
  colnames(Metrics) <- paste0(condition_name,"_",colnames(Metrics))
  Mean_Var_Table <- cbind(Mean_Var_Table,Metrics)
}

View(VST_Metric_Correlation_Table)
View(Mean_Var_Table)
# Very high mean-variability correlation coefficients - may as well just use VST() MADs


# Plot these out to demonstrate that the pattern is similar...
# Variability metrics against log CPM
# Adj Var ~ logCPM
Ctrl_MeanVar_AdjVar_logCPM_plots <- Variability_mean_plots_precomp(MeanVarTable = Mean_Var_Table,log_mean = F,log_var = T,mean_name = "Mean logCPM (Ctrl)",var_name = "log10(Adjusted variance)",mean_col = 'Ctrl_logCPM_Mean',var_col ='Ctrl_Adj_Var')
Ctrl_MeanVar_AdjVar_logCPM_plots_full <- Ctrl_MeanVar_AdjVar_logCPM_plots[[1]]
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = Ctrl_MeanVar_AdjVar_logCPM_plots_full,filename = 'Ctrl_MeanVar_AdjVar_logCPM_plots_full.svg',dpi=300,width=2.24,height = 2.24)
HS_MeanVar_AdjVar_logCPM_plots <- Variability_mean_plots_precomp(MeanVarTable = Mean_Var_Table,log_mean = F,log_var = T,mean_name = "Mean logCPM (HS)",var_name = "log10(Adjusted variance)",mean_col = 'HS_logCPM_Mean',var_col ='HS_Adj_Var')
HS_MeanVar_AdjVar_logCPM_plots_full <- HS_MeanVar_AdjVar_logCPM_plots[[1]]
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = HS_MeanVar_AdjVar_logCPM_plots_full,filename = 'HS_MeanVar_AdjVar_logCPM_plots_full.svg',dpi=300,width=2.24,height = 2.24)

# LOESS Var ~ logCPM
Ctrl_MeanVar_LOESSVar_logCPM_plots <- Variability_mean_plots_precomp(MeanVarTable = Mean_Var_Table,log_mean = F,log_var = F,mean_name = "Mean logCPM (Ctrl)",var_name = "LOESS variance",mean_col = 'Ctrl_logCPM_Mean',var_col ='Ctrl_LOESS_Var')
Ctrl_MeanVar_LOESSVar_logCPM_plots_full <- Ctrl_MeanVar_LOESSVar_logCPM_plots[[1]]
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = Ctrl_MeanVar_LOESSVar_logCPM_plots_full,filename = 'Ctrl_MeanVar_LOESSVar_logCPM_plots_full.svg',dpi=300,width=2.24,height = 2.24)
HS_MeanVar_LOESSVar_logCPM_plots <- Variability_mean_plots_precomp(MeanVarTable = Mean_Var_Table,log_mean = F,log_var = F,mean_name = "Mean logCPM (HS)",var_name = "LOESS variance",mean_col = 'HS_logCPM_Mean',var_col ='HS_LOESS_Var')
HS_MeanVar_LOESSVar_logCPM_plots_full <- HS_MeanVar_LOESSVar_logCPM_plots[[1]]
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = HS_MeanVar_LOESSVar_logCPM_plots_full,filename = 'HS_MeanVar_LOESSVar_logCPM_plots_full.svg',dpi=300,width=2.24,height = 2.24)


####################################
######### End) Save files  #########
####################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
ggsave(plot = QC_plots,filename = "QC_plots_RNAseq.png",height=6,width=6)
write.csv(VST_Metric_Correlation_Table,"Different_VST_Metrics_Correlation_Table.csv",row.names = F)
write.csv(Mean_Var_Table,"Gene_MeanVar_Table.csv",row.names = F)

save.image(file='Decide_MeanVar_Metrics.RData')
