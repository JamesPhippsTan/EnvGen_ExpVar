# Investigating eQTL and veQTL 
# Part 3: Null distribution generation for derived allele frequencies of different categories
# Gets the median DAF and skewness of each distribution 

# Last Updated: 11/9/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls(all = T))

library(dplyr)
library(ggplot2)
library(moments)

# Continues with objects saved in the following .RData files
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
load(file='5_Investigating_eQTL_and_veQTL.RData')
load(file='6_Investigating_eQTL_and_veQTL_Allele_Age1.RData')
load(file='6_Investigating_eQTL_and_veQTL_Allele_Age2.RData')

######################################################
# (6) Figuring out the random expectation for DAFs ###
######################################################

# Method 1 - using the MAF information but not the increased or decreased information 

# Save files to here
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\DAF_downsampling")

# HS trans-v 
# derived increasing 
HS_trans_veQTL_der_inc_MAF <- MAF_DAF[HS_trans_veQTL_der_inc$SNP,'HS_MAF']
HS_trans_veQTL_der_inc_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  HS_trans_veQTL_der_inc_MAF_Criteria[quintile,'Min']=quantile(HS_trans_veQTL_der_inc_MAF,min)
  HS_trans_veQTL_der_inc_MAF_Criteria[quintile,'Max']=quantile(HS_trans_veQTL_der_inc_MAF,max)
  HS_trans_veQTL_der_inc_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
View(HS_trans_veQTL_der_inc_MAF_Criteria)
# derived decreasing 
HS_trans_veQTL_der_dec_MAF <- MAF_DAF[HS_trans_veQTL_der_dec$SNP,'HS_MAF']
HS_trans_veQTL_der_dec_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  HS_trans_veQTL_der_dec_MAF_Criteria[quintile,'Min']=quantile(HS_trans_veQTL_der_dec_MAF,min)
  HS_trans_veQTL_der_dec_MAF_Criteria[quintile,'Max']=quantile(HS_trans_veQTL_der_dec_MAF,max)
  HS_trans_veQTL_der_dec_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
# What dataframes are necessary
HS_tv_DI_MAF_Criteria <- HS_trans_veQTL_der_inc_MAF_Criteria # Criteria 1 
HS_tv_DD_MAF_Criteria <- HS_trans_veQTL_der_dec_MAF_Criteria # Criteria 2
HS_tv_DI_df <- make_SNP_MAF_DAF_df(MAF_DAF[HS_trans_veQTL_der_inc$SNP,c('HS_MAF','HS_DAF')]) # Subsample der_inc
HS_tv_DD_df <- make_SNP_MAF_DAF_df(MAF_DAF[HS_trans_veQTL_der_dec$SNP,c('HS_MAF','HS_DAF')]) # Subsample der_dec
HS_tv_DI_non_df <- make_SNP_MAF_DAF_df(MAF_DAF[non_HS_veQTL_aged,c('HS_MAF','HS_DAF')]) # Subsample non-veQTL
HS_tv_DD_non_df <- make_SNP_MAF_DAF_df(MAF_DAF[non_HS_veQTL_aged,c('HS_MAF','HS_DAF')]) # Subsample non-veQTL
# Sufficient information for downsampling and estimation of average DAF...

# HS trans-e
# deried increasing 
HS_trans_eQTL_der_inc_MAF <- MAF_DAF[HS_trans_eQTL_der_inc$SNP,'HS_MAF']
HS_trans_eQTL_der_inc_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  HS_trans_eQTL_der_inc_MAF_Criteria[quintile,'Min']=quantile(HS_trans_eQTL_der_inc_MAF,min)
  HS_trans_eQTL_der_inc_MAF_Criteria[quintile,'Max']=quantile(HS_trans_eQTL_der_inc_MAF,max)
  HS_trans_eQTL_der_inc_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
View(HS_trans_eQTL_der_inc_MAF_Criteria)
# deried decreasing 
HS_trans_eQTL_der_dec_MAF <- MAF_DAF[HS_trans_eQTL_der_dec$SNP,'HS_MAF']
HS_trans_eQTL_der_dec_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  HS_trans_eQTL_der_dec_MAF_Criteria[quintile,'Min']=quantile(HS_trans_eQTL_der_dec_MAF,min)
  HS_trans_eQTL_der_dec_MAF_Criteria[quintile,'Max']=quantile(HS_trans_eQTL_der_dec_MAF,max)
  HS_trans_eQTL_der_dec_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
# What dataframes are necessary
HS_te_DI_MAF_Criteria <- HS_trans_eQTL_der_inc_MAF_Criteria # Criteria 1 
HS_te_DD_MAF_Criteria <- HS_trans_eQTL_der_dec_MAF_Criteria # Criteria 2
HS_te_DI_df <- make_SNP_MAF_DAF_df(MAF_DAF[HS_trans_eQTL_der_inc$SNP,c('HS_MAF','HS_DAF')]) # Subsample der_inc
HS_te_DD_df <- make_SNP_MAF_DAF_df(MAF_DAF[HS_trans_eQTL_der_dec$SNP,c('HS_MAF','HS_DAF')]) # Subsample der_dec
HS_te_DI_non_df <- make_SNP_MAF_DAF_df(MAF_DAF[non_HS_trans_eQTL_by_MAF$SNP,c('HS_MAF','HS_DAF')]) # Subsample non-eQTL
HS_te_DD_non_df <- make_SNP_MAF_DAF_df(MAF_DAF[non_HS_trans_eQTL_by_MAF$SNP,c('HS_MAF','HS_DAF')]) # Subsample non-eQTL
# Sufficient information for downsampling and estimation of aerage DAF...

# Ctrl trans-v 
# derived increasing 
Ctrl_trans_veQTL_der_inc_MAF <- MAF_DAF[Ctrl_trans_veQTL_der_inc$SNP,'Ctrl_MAF']
Ctrl_trans_veQTL_der_inc_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  Ctrl_trans_veQTL_der_inc_MAF_Criteria[quintile,'Min']=quantile(Ctrl_trans_veQTL_der_inc_MAF,min)
  Ctrl_trans_veQTL_der_inc_MAF_Criteria[quintile,'Max']=quantile(Ctrl_trans_veQTL_der_inc_MAF,max)
  Ctrl_trans_veQTL_der_inc_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
View(Ctrl_trans_veQTL_der_inc_MAF_Criteria)
# derived decreasing 
Ctrl_trans_veQTL_der_dec_MAF <- MAF_DAF[Ctrl_trans_veQTL_der_dec$SNP,'Ctrl_MAF']
Ctrl_trans_veQTL_der_dec_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  Ctrl_trans_veQTL_der_dec_MAF_Criteria[quintile,'Min']=quantile(Ctrl_trans_veQTL_der_dec_MAF,min)
  Ctrl_trans_veQTL_der_dec_MAF_Criteria[quintile,'Max']=quantile(Ctrl_trans_veQTL_der_dec_MAF,max)
  Ctrl_trans_veQTL_der_dec_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
# What dataframes are necessary
Ctrl_tv_DI_MAF_Criteria <- Ctrl_trans_veQTL_der_inc_MAF_Criteria # Criteria 1 
Ctrl_tv_DD_MAF_Criteria <- Ctrl_trans_veQTL_der_dec_MAF_Criteria # Criteria 2
Ctrl_tv_DI_df <- make_SNP_MAF_DAF_df(MAF_DAF[Ctrl_trans_veQTL_der_inc$SNP,c('Ctrl_MAF','Ctrl_DAF')]) # Subsample der_inc
Ctrl_tv_DD_df <- make_SNP_MAF_DAF_df(MAF_DAF[Ctrl_trans_veQTL_der_dec$SNP,c('Ctrl_MAF','Ctrl_DAF')]) # Subsample der_dec
Ctrl_tv_DI_non_df <- make_SNP_MAF_DAF_df(MAF_DAF[non_Ctrl_veQTL_aged,c('Ctrl_MAF','Ctrl_DAF')]) # Subsample non-veQTL
Ctrl_tv_DD_non_df <- make_SNP_MAF_DAF_df(MAF_DAF[non_Ctrl_veQTL_aged,c('Ctrl_MAF','Ctrl_DAF')]) # Subsample non-veQTL
# Sufficient information for downsampling and estimation of average DAF...
# Remove sampling from table step
# Sample rows 1000 times and compute average DAF per sample
# SNPs chosen are non-overlapping

# Ctrl trans-e
# derived increasing 
Ctrl_trans_eQTL_der_inc_MAF <- MAF_DAF[Ctrl_trans_eQTL_der_inc$SNP,'Ctrl_MAF']
Ctrl_trans_eQTL_der_inc_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  Ctrl_trans_eQTL_der_inc_MAF_Criteria[quintile,'Min']=quantile(Ctrl_trans_eQTL_der_inc_MAF,min)
  Ctrl_trans_eQTL_der_inc_MAF_Criteria[quintile,'Max']=quantile(Ctrl_trans_eQTL_der_inc_MAF,max)
  Ctrl_trans_eQTL_der_inc_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
View(Ctrl_trans_eQTL_der_inc_MAF_Criteria)
# deried decreasing 
Ctrl_trans_eQTL_der_dec_MAF <- MAF_DAF[Ctrl_trans_eQTL_der_dec$SNP,'Ctrl_MAF']
Ctrl_trans_eQTL_der_dec_MAF_Criteria <- data.frame(Quintile=1:5,Min=1:5,Max=1:5,NumSNPstoSample=1:5)
for (quintile in 1:5){
  min = (quintile-1)*0.2
  max = 0.2+min
  Ctrl_trans_eQTL_der_dec_MAF_Criteria[quintile,'Min']=quantile(Ctrl_trans_eQTL_der_dec_MAF,min)
  Ctrl_trans_eQTL_der_dec_MAF_Criteria[quintile,'Max']=quantile(Ctrl_trans_eQTL_der_dec_MAF,max)
  Ctrl_trans_eQTL_der_dec_MAF_Criteria[quintile,'NumSNPstoSample']=1000
}
# What dataframes are necessary
Ctrl_te_DI_MAF_Criteria <- Ctrl_trans_eQTL_der_inc_MAF_Criteria # Criteria 1 
Ctrl_te_DD_MAF_Criteria <- Ctrl_trans_eQTL_der_dec_MAF_Criteria # Criteria 2
Ctrl_te_DI_df <- make_SNP_MAF_DAF_df(MAF_DAF[Ctrl_trans_eQTL_der_inc$SNP,c('Ctrl_MAF','Ctrl_DAF')]) # Subsample der_inc
Ctrl_te_DD_df <- make_SNP_MAF_DAF_df(MAF_DAF[Ctrl_trans_eQTL_der_dec$SNP,c('Ctrl_MAF','Ctrl_DAF')]) # Subsample der_dec
Ctrl_te_DI_non_df <- make_SNP_MAF_DAF_df(MAF_DAF[non_Ctrl_trans_eQTL_by_MAF$SNP,c('Ctrl_MAF','Ctrl_DAF')]) # Subsample non-eQTL
Ctrl_te_DD_non_df <- make_SNP_MAF_DAF_df(MAF_DAF[non_Ctrl_trans_eQTL_by_MAF$SNP,c('Ctrl_MAF','Ctrl_DAF')]) # Subsample non-eQTL
# Sufficient information for downsampling and estimation of average DAF...
# Remove sampling from table step
# Sample rows 1000 times and compute aerage DAF per sample
# SNPs chosen are non-overlapping

# Save the files
# Ctrl trans-veQTL
write.table(Ctrl_tv_DI_MAF_Criteria,'Ctrl_tv_DI_MAF_Criteria.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_tv_DD_MAF_Criteria,'Ctrl_tv_DD_MAF_Criteria.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_tv_DI_df,'Ctrl_tv_DI_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_tv_DD_df,'Ctrl_tv_DD_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_tv_DI_non_df,'Ctrl_tv_DI_non_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_tv_DD_non_df,'Ctrl_tv_DD_non_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
# HS trans-veQTL
write.table(HS_tv_DI_MAF_Criteria,'HS_tv_DI_MAF_Criteria.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_tv_DD_MAF_Criteria,'HS_tv_DD_MAF_Criteria.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_tv_DI_df,'HS_tv_DI_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_tv_DD_df,'HS_tv_DD_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_tv_DI_non_df,'HS_tv_DI_non_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_tv_DD_non_df,'HS_tv_DD_non_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
# Ctrl trans-eQTL
write.table(Ctrl_te_DI_MAF_Criteria,'Ctrl_te_DI_MAF_Criteria.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_te_DD_MAF_Criteria,'Ctrl_te_DD_MAF_Criteria.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_te_DI_df,'Ctrl_te_DI_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_te_DD_df,'Ctrl_te_DD_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_te_DI_non_df,'Ctrl_te_DI_non_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(Ctrl_te_DD_non_df,'Ctrl_te_DD_non_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
# HS trans-eQTL
write.table(HS_te_DI_MAF_Criteria,'HS_te_DI_MAF_Criteria.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_te_DD_MAF_Criteria,'HS_te_DD_MAF_Criteria.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_te_DI_df,'HS_te_DI_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_te_DD_df,'HS_te_DD_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_te_DI_non_df,'HS_te_DI_non_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 
write.table(HS_te_DD_non_df,'HS_te_DD_non_downsubsample_SNPs.txt',row.names = F,col.names = T,quote = F) 


###################################################################
# (7) Getting summary statistics of the downsampled distributions #
###################################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\DAF_downsampling")

#------------------------------------------------------------
# Helper: subsample SNPs given MAF criteria
#------------------------------------------------------------
subsample_snps <- function(snp_df, criteria) {
  do.call(rbind, lapply(1:nrow(criteria), function(i) {
    crit <- criteria[i, ]
    eligible <- snp_df %>% filter(MAF >= crit$Min, MAF < crit$Max)
    if (nrow(eligible) < crit$NumSNPstoSample) {
      stop(paste("Not enough SNPs in range for quintile", crit$Quintile))
    }
    eligible %>% slice_sample(n = crit$NumSNPstoSample)
  }))
}

#------------------------------------------------------------
# Compute one replicate: skew, median, KS p-value (shared subsampling)
#------------------------------------------------------------
compute_skew_median_once <- function(qtl_df, nonqtl_df, criteria, ctrl_fixed_qtl = FALSE, qtl_fixed_daf = NULL) {
  # If ctrl_fixed_qtl is TRUE, qtl_df won't be subsampled; qtl_fixed_daf must be provided (vector)
  if (ctrl_fixed_qtl) {
    qtl_sub_daf <- qtl_fixed_daf
  } else {
    qtl_sub <- subsample_snps(qtl_df, criteria)
    qtl_sub_daf <- qtl_sub$DAF
  }
  nonqtl_sub <- subsample_snps(nonqtl_df, criteria)
  nonqtl_sub_daf <- nonqtl_sub$DAF
  
  # compute skew and median
  skew_qtl <- skewness(qtl_sub_daf, na.rm = TRUE)
  skew_nonqtl <- skewness(nonqtl_sub_daf, na.rm = TRUE)
  median_qtl <- median(qtl_sub_daf, na.rm = TRUE)
  median_nonqtl <- median(nonqtl_sub_daf, na.rm = TRUE)
  
  # KS test (safe with tryCatch)
  ks_p <- tryCatch(ks.test(qtl_sub_daf, nonqtl_sub_daf)$p.value, error = function(e) NA)
  
  data.frame(
    skew_qtl = skew_qtl,
    skew_nonqtl = skew_nonqtl,
    median_qtl = median_qtl,
    median_nonqtl = median_nonqtl,
    ks_pvalue = ks_p
  )
}

#------------------------------------------------------------
# Replicate wrapper (shared subsampling)
#------------------------------------------------------------
sample_skew_median_replicates <- function(qtl_df, nonqtl_df, criteria, nrep = 1000, ctrl_fixed_qtl = FALSE, qtl_fixed_daf = NULL) {
  bind_rows(lapply(1:nrep, function(r) {
    set.seed(r)
    compute_skew_median_once(qtl_df, nonqtl_df, criteria, ctrl_fixed_qtl, qtl_fixed_daf)
  }))
}

#------------------------------------------------------------
# Plot combined DAF distributions (3 reps)
#------------------------------------------------------------
plot_combined_daf <- function(qtl_df, nonqtl_df, criteria, cat, nrep = 3, ctrl_fixed_qtl = FALSE, qtl_fixed_daf = NULL) {
  for (r in 1:nrep) {
    set.seed(r)
    if (ctrl_fixed_qtl) {
      # use the fixed qtl_daf vector as-is; wrap into a data.frame for plotting
      qtl_sub <- data.frame(DAF = qtl_fixed_daf, Type = "QTL")
    } else {
      qtl_sub <- subsample_snps(qtl_df, criteria) %>% mutate(Type = "QTL")
    }
    nonqtl_sub <- subsample_snps(nonqtl_df, criteria) %>% mutate(Type = "nonQTL")
    subsample_all <- bind_rows(qtl_sub, nonqtl_sub)
    
    p_combined <- ggplot(subsample_all, aes(x = DAF, fill = Type)) +
      geom_histogram(position = "identity", alpha = 0.5, bins = 10) +
      theme_classic() +
      scale_fill_manual(values = c("QTL" = "skyblue", "nonQTL" = "orange")) +
      labs(x = "DAF", y = "Count", title = paste("QTL vs non-QTL - Rep", r))
    
    ggsave(filename = paste0(cat, "_rep", r, "_Combined_DAF_Distribution.svg"),
           plot = p_combined, width = 4, height = 3, dpi = 300)
  }
}

#------------------------------------------------------------
# Plot non-QTL only (used when QTL fixed / cannot be subsampled)
#------------------------------------------------------------
plot_nonqtl_daf <- function(snp_df, criteria, cat, nrep = 3) {
  for (r in 1:nrep) {
    set.seed(r)
    subsample <- subsample_snps(snp_df, criteria)
    p <- ggplot(subsample, aes(x = DAF)) +
      geom_histogram(bins = 10, fill = "skyblue") +
      theme_classic() +
      labs(x = "DAF", y = "Count", title = paste("non-QTL SNPs - Rep", r))
    
    ggsave(filename = paste0(cat, "_rep", r, "_nonQTL_DAF_Distribution.svg"),
           plot = p, width = 3, height = 3, dpi = 300)
  }
}

#------------------------------------------------------------
# Main loop
#------------------------------------------------------------
categories <- c("HS_te_DI", "HS_te_DD", 
                "Ctrl_te_DI", "Ctrl_te_DD",
                "HS_tv_DI", "HS_tv_DD",
                "Ctrl_tv_DI", "Ctrl_tv_DD")

summary_table <- data.frame(
  Category = character(),
  mean_QTL_skew = numeric(),
  mean_nonQTL_skew = numeric(),
  skew_t_pvalue = numeric(),
  mean_QTL_medianDAF = numeric(),
  mean_nonQTL_medianDAF = numeric(),
  median_t_pvalue = numeric(),
  mean_KS_pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (cat in categories) {
  message("Processing ", cat, " ...")
  criteria <- read.table(paste0(cat, "_MAF_Criteria.txt"), header = TRUE)
  qtl      <- read.table(paste0(cat, "_downsubsample_SNPs.txt"), header = TRUE)
  nonqtl   <- read.table(paste0(cat, "_non_downsubsample_SNPs.txt"), header = TRUE)
  set.seed(4)
  
  # handle Ctrl_tv cases where QTL cannot be subsampled (use actual QTL DAFs)
  if (cat %in% c("Ctrl_tv_DI", "Ctrl_tv_DD")) {
    qtl_fixed_daf <- qtl$DAF
    # Plot nonQTL distributions (3 reps)
    plot_nonqtl_daf(nonqtl, criteria, cat, nrep = 3)
    # For summary & replicate results: compute 1000 replicates where QTL is fixed and nonQTL subsampled
    df_reps <- bind_rows(lapply(1:1000, function(r) {
      set.seed(r)
      non_sub <- subsample_snps(nonqtl, criteria)
      ks_p <- tryCatch(ks.test(qtl_fixed_daf, non_sub$DAF)$p.value, error = function(e) NA)
      data.frame(
        skew_qtl = skewness(qtl_fixed_daf, na.rm = TRUE),
        skew_nonqtl = skewness(non_sub$DAF, na.rm = TRUE),
        median_qtl = median(qtl_fixed_daf, na.rm = TRUE),
        median_nonqtl = median(non_sub$DAF, na.rm = TRUE),
        ks_pvalue = ks_p
      )
    }))
    df_reps$replicate <- 1:1000
    
    # t-tests between the replicated vectors (QTL repeated vs nonQTL reps)
    t_res_skew <- tryCatch(wilcox.test(df_reps$skew_qtl, df_reps$skew_nonqtl), error = function(e) NULL)
    t_res_median <- tryCatch(wilcox.test(df_reps$median_qtl, df_reps$median_nonqtl), error = function(e) NULL)
  } else {
    # normal case: shared subsampling inside replicate function -> returns per-replicate ks p-values too
    df_reps <- sample_skew_median_replicates(qtl, nonqtl, criteria, nrep = 1000, ctrl_fixed_qtl = FALSE, qtl_fixed_daf = NULL)
    df_reps$replicate <- 1:1000
    # plot 3 replicates
    plot_combined_daf(qtl, nonqtl, criteria, cat, nrep = 3, ctrl_fixed_qtl = FALSE)
    # tests
    t_res_skew <- tryCatch(wilcox.test(df_reps$skew_qtl, df_reps$skew_nonqtl), error = function(e) NULL)
    t_res_median <- tryCatch(wilcox.test(df_reps$median_qtl, df_reps$median_nonqtl), error = function(e) NULL)
  }
  
  # Save replicate-level results (skew + median + ks)
  out_rep_file <- paste0(cat, "_skew_median_results.txt")
  write.table(df_reps, out_rep_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Build summary values
  mean_qtl_skew <- mean(df_reps$skew_qtl, na.rm = TRUE)
  mean_nonqtl_skew <- mean(df_reps$skew_nonqtl, na.rm = TRUE)
  mean_qtl_med <- mean(df_reps$median_qtl, na.rm = TRUE)
  mean_nonqtl_med <- mean(df_reps$median_nonqtl, na.rm = TRUE)
  mean_ks <- mean(df_reps$ks_pvalue, na.rm = TRUE)
  
  summary_table <- rbind(summary_table, data.frame(
    Category = cat,
    mean_QTL_skew = mean_qtl_skew,
    mean_nonQTL_skew = mean_nonqtl_skew,
    skew_t_pvalue = ifelse(is.null(t_res_skew), NA, t_res_skew$p.value),
    mean_QTL_medianDAF = mean_qtl_med,
    mean_nonQTL_medianDAF = mean_nonqtl_med,
    median_t_pvalue = ifelse(is.null(t_res_median), NA, t_res_median$p.value),
    mean_KS_pvalue = mean_ks,
    stringsAsFactors = FALSE
  ))
  
  message("Saved replicate file: ", out_rep_file)
}

# Save combined summary table
write.table(summary_table, "DAF_Downsampling_Combined_Summary_Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
