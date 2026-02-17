# Investigating eQTL and veQTL - only the allele age component so I can run as a background job.
# Part 1: Just the per-QTL FDI

# Last Updated: 17/9/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls(all = T))

library(reshape)
library(tidyr)
library(dplyr)
library(stats)
library(EnvStats)
library(data.table)
library(ggplot2)
library(UpSetR)
library(tibble)
library(ggvenn)
library(gridExtra)
library(patchwork)
library(eulerr)
library(svglite)

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
 load(file='6_Investigating_eQTL_and_veQTL_Allele_Age1.RData')

# Functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('gProfiler_Functions.R')
source('Quantile_Functions.R')
source('Variability_Functions.R')
source('QTL_Analysis_Functions.R') 

#######################
##### Datasets ########
#######################

# Continues with objects saved in the following .RData file
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
load(file='5_Investigating_eQTL_and_veQTL.RData')

#################################################################################################################
# (1) Are older or newer alleles more likely to increase transcript level variability? eQTL and veQTL.... #######
#################################################################################################################

# Define whether variability increased or decreased for each derived vSNP allele. 
# Repeat the same thing for multiple sets of non-significant SNPs.
# Then run a comparison of the Fractions of the genes regulated by a Derived vSNPs that Increase 
# FDI for short

# What are the names of the columns harbouring the direction of change? 
# ciseqtl=slope,transeqtl=b,cisveqtl=COR,transveqtl=COR
# What are the names of the columns harbouring the SNPs? Which one to prioritize?
# A1 and A2 - all dataframes; REF and ALT are also present veQTL mapping


#########
#  eQTL #
#########

Ctrl_cis_per_eQTL_FDI <- get_SNP_specific_FDI(QTL_df = Ctrl_cis_eqtl_sig,
                                              SNP_ages_df = SNP_allele_age,
                                              slope_allele_column = 'A2',
                                              slope_column = 'slope')

Ctrl_trans_per_eQTL_FDI <- get_SNP_specific_FDI(QTL_df = Ctrl_trans_eqtl_sig,
                                                SNP_ages_df = SNP_allele_age,
                                                slope_allele_column = 'A2',
                                                slope_column = 'b')

HS_cis_per_eQTL_FDI <- get_SNP_specific_FDI(QTL_df = HS_cis_eqtl_sig,
                                            SNP_ages_df = SNP_allele_age,
                                            slope_allele_column = 'A2',
                                            slope_column = 'slope')

HS_trans_per_eQTL_FDI <- get_SNP_specific_FDI(QTL_df = HS_trans_eqtl_sig,
                                              SNP_ages_df = SNP_allele_age,
                                              slope_allele_column = 'A2',
                                              slope_column = 'b')


#########
# veQTL #
#########

# This is because some SNPs can regulate many genes

Ctrl_cis_per_veQTL_FDI <- get_SNP_specific_FDI(QTL_df = Ctrl_cis_veQTL_sig,
                                              SNP_ages_df = SNP_allele_age,
                                              slope_allele_column = 'ALT',
                                              slope_column = 'COR')

Ctrl_trans_per_veQTL_FDI <- get_SNP_specific_FDI(QTL_df = Ctrl_trans_veQTL_sig,
                                                SNP_ages_df = SNP_allele_age,
                                                slope_allele_column = 'ALT',
                                                slope_column = 'COR')

HS_cis_per_veQTL_FDI <- get_SNP_specific_FDI(QTL_df = HS_cis_veQTL_sig,
                                            SNP_ages_df = SNP_allele_age,
                                            slope_allele_column = 'ALT',
                                            slope_column = 'COR')

HS_trans_per_veQTL_FDI <- get_SNP_specific_FDI(QTL_df = HS_trans_veQTL_sig,
                                              SNP_ages_df = SNP_allele_age,
                                              slope_allele_column = 'ALT',
                                              slope_column = 'COR')

# Save tables for FDI per SNP at 100% and at 0%
# List all your dataframes (replace with your actual dataframe names)
# Method 1: Using a list of dataframes (recommended)
df_list <- list(Ctrl_cis_eqtl=Ctrl_cis_per_eQTL_FDI, 
                HS_cis_eqtl=HS_cis_per_eQTL_FDI, 
                Ctrl_trans_eqtl=Ctrl_trans_per_eQTL_FDI, 
                HS_trans_eqtl=HS_trans_per_eQTL_FDI,
                Ctrl_cis_veqtl=Ctrl_cis_per_veQTL_FDI, 
                HS_cis_veqtl=HS_cis_per_veQTL_FDI, 
                Ctrl_trans_veqtl=Ctrl_trans_per_veQTL_FDI, 
                HS_trans_veqtl=HS_trans_per_veQTL_FDI)

summary_df <- data.frame(
  row.names = names(df_list),
  Number_of_SNPs = sapply(df_list, function(df) nrow(df)),
  Increase_only_count = sapply(df_list, function(df) sum(df$FDI == 1, na.rm = TRUE)),  # Raw count (FDI=1)
  Increase_only_percent = sapply(df_list, function(df) mean(df$FDI == 1, na.rm = TRUE) * 100),  # Percentage (FDI=1)
  Decrease_only_count = sapply(df_list, function(df) sum(df$FDI == 0, na.rm = TRUE)),  # Raw count (FDI=0)
  Decrease_only_percent = sapply(df_list, function(df) mean(df$FDI == 0, na.rm = TRUE) * 100)  # Percentage (FDI=0)
)
summary_df$Both_count <- summary_df$Number_of_SNPs-summary_df$Increase_only_count-summary_df$Decrease_only_count
summary_df$Both_percent <- 100-summary_df$Increase_only_percent-summary_df$Decrease_only_percent
summary_df$Increase_only_percent_noboth <- 100*summary_df$Increase_only_count/(summary_df$Increase_only_count+summary_df$Decrease_only_count)
summary_df$Decrease_only_percent_noboth <- 100*summary_df$Decrease_only_count/(summary_df$Increase_only_count+summary_df$Decrease_only_count)

# View the result
summary_df[,2:9]<- round(summary_df[,2:9],digits=1)
rownames(summary_df) <- gsub('_',' ',rownames(summary_df))
rownames(summary_df) <- gsub('cis ','cis-',rownames(summary_df))
rownames(summary_df) <- gsub('trans ','trans-',rownames(summary_df))
View(summary_df)

# Save the result
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Derived_Allele_Increased_Fraction")
write.csv(summary_df,file = 'FDI_FDD_percent_base.csv',row.names = T)

############################
# Plotting the percentages #
############################

# Check the necessary numbers
rownames(summary_df) <- names(df_list)

# Horizontal barplot visualisation
data_trans_eQTL_Ctrl <- data.frame(
  Partition = c("Increases only", "Decreases only", "Both"),
  Value = c(summary_df['Ctrl_trans_eqtl','Increase_only_count'], 
            summary_df['Ctrl_trans_eqtl','Decrease_only_count'],
            summary_df['Ctrl_trans_eqtl','Both_count']) 
) %>%
  mutate(
    Percent = Value / sum(Value),
    Partition = factor(Partition, levels = c("Both","Decreases only","Increases only"))
  )

data_trans_eQTL_HS <- data.frame(
  Partition = c("Increases only", "Decreases only", "Both"),
  Value = c(summary_df['HS_trans_eqtl','Increase_only_count'], 
            summary_df['HS_trans_eqtl','Decrease_only_count'],
            summary_df['HS_trans_eqtl','Both_count']) 
) %>%
  mutate(
    Percent = Value / sum(Value),
    Partition = factor(Partition, levels = c("Both","Decreases only","Increases only"))
  )

data_trans_veQTL_Ctrl <- data.frame(
  Partition = c("Increases only", "Decreases only", "Both"),
  Value = c(summary_df['Ctrl_trans_veqtl','Increase_only_count'], 
            summary_df['Ctrl_trans_veqtl','Decrease_only_count'],
            summary_df['Ctrl_trans_veqtl','Both_count']) 
) %>%
  mutate(
    Percent = Value / sum(Value),
    Partition = factor(Partition, levels = c("Both","Decreases only","Increases only"))
  )

data_trans_veQTL_HS <- data.frame(
  Partition = c("Increases only", "Decreases only", "Both"),
  Value = c(summary_df['HS_trans_veqtl','Increase_only_count'], 
            summary_df['HS_trans_veqtl','Decrease_only_count'],
            summary_df['HS_trans_veqtl','Both_count']) 
) %>%
  mutate(
    Percent = Value / sum(Value),
    Partition = factor(Partition, levels = c("Both","Decreases only","Increases only"))
  )


# Colours
eQTL_colors <- c(
  "Both" = "#4F8E4D",
  "Decreases only" = "#77BB75",
  "Increases only" = "#2B5D29"
)
veQTL_colors <- c(
  "Both" = "#9651EC",
  "Decreases only" = "#CFB2F3",
  "Increases only" = "#611BB8"
)

# Create overlap plots
trans_eQTL_Ctrl_derplot <- create_overlap_plot(data_trans_eQTL_Ctrl,eQTL_colors,
                                     '# derived eQTL alleles (Ctrl)')
trans_veQTL_Ctrl_derplot <- create_overlap_plot(data_trans_veQTL_Ctrl,veQTL_colors,
                                               '# derived veQTL alleles (Ctrl)')
trans_eQTL_HS_derplot <- create_overlap_plot(data_trans_eQTL_HS,eQTL_colors,
                                               '# derived eQTL alleles (HS)')
trans_veQTL_HS_derplot <- create_overlap_plot(data_trans_veQTL_HS,veQTL_colors,
                                               '# derived veQTL alleles (HS)')


setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Derived_Allele_Increased_Fraction")
ggsave(plot = trans_eQTL_Ctrl_derplot, filename='trans_eQTL_Ctrl_derplot.svg', width = 4, height = 1.5, dpi = 300)
ggsave(plot = trans_veQTL_Ctrl_derplot, filename='trans_veQTL_Ctrl_derplot.svg', width = 4, height =1.5, dpi = 300)
ggsave(plot = trans_eQTL_HS_derplot, filename='trans_eQTL_HS_derplot.svg', width = 4, height = 1.5, dpi = 300)
ggsave(plot = trans_veQTL_HS_derplot, filename='trans_veQTL_HS_derplot.svg', width = 4, height =1.5, dpi = 300)


#################################################
##### End) Save working environment #############
#################################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
save.image(file='6_Investigating_eQTL_and_veQTL_Allele_Age1.RData')

# Old visualisation - histograms #

Ctrl_cis_veQTL_FDI_plot <- ggplot(Ctrl_cis_per_veQTL_FDI, aes(x = 100*(1-FDI))) + 
  annotate("text", y = Inf, x = 50, label = paste0("Total number of veQTLs = ",nrow(Ctrl_cis_per_veQTL_FDI)), size = 4, hjust = 0.5,vjust = 1)+
  geom_histogram(bins=5, fill = "lightblue", color = "black", alphsig = 0.7) + 
  labs(x = "Percent of Derived Alleles Decreasing Variability", 
       y = "Number of cis-veQTLs - Ctrl") + theme_classic()
Ctrl_cis_veQTL_FDI_plot

Ctrl_trans_veQTL_FDI_plot <- ggplot(Ctrl_trans_per_veQTL_FDI, aes(x = 100*(1-FDI))) + 
  annotate("text", y = Inf, x = 50, label = paste0("Total number of veQTLs = ",nrow(Ctrl_trans_per_veQTL_FDI)), size = 4, hjust = 0.5,vjust = 1)+
  geom_histogram(bins=5, fill = "lightblue", color = "black", alphsig = 0.7) + 
  labs(x = "Percent of Derived Alleles Decreasing Variability", 
       y = "Number of trans-veQTLs - Ctrl") + theme_classic()
Ctrl_trans_veQTL_FDI_plot

HS_cis_veQTL_FDI_plot <- ggplot(HS_cis_per_veQTL_FDI, aes(x = 100*(1-FDI))) + 
  annotate("text", y = Inf, x = 100, label = paste0("Total number of veQTLs = ",nrow(HS_cis_per_veQTL_FDI)), size = 4, hjust = 0.5,vjust = 1)+
  geom_histogram(bins=5, fill = "lightblue", color = "black", alphsig = 0.7) + 
  labs(x = "Percent of Derived Alleles Decreasing Variability", 
       y = "Number of cis-veQTLs - HS") + theme_classic()
HS_cis_veQTL_FDI_plot

HS_trans_veQTL_FDI_plot <- ggplot(HS_trans_per_veQTL_FDI, aes(x = 100*(1-FDI))) + 
  annotate("text", y = Inf, x = 50, label = paste0("Total number of veQTLs = ",nrow(HS_trans_per_veQTL_FDI)), size = 4, hjust = 0.5,vjust = 1)+
  geom_histogram(bins=5, fill = "lightblue", color = "black", alphsig = 0.7) + 
  labs(x = "Percent of Derived Alleles Decreasing Variability", 
       y = "Number of trans-veQTLs - HS") + theme_classic()
HS_trans_veQTL_FDI_plot


Ctrl_cis_eQTL_FDI_plot <- ggplot(Ctrl_cis_per_eQTL_FDI, aes(x = 100*(1-FDI))) + 
  annotate("text", y = Inf, x = 50, label = paste0("Total number of eQTLs = ",nrow(Ctrl_cis_per_eQTL_FDI)), size = 4, hjust = 0.5,vjust = 1)+
  geom_histogram(bins=5, fill = "lightblue", color = "black", alphsig = 0.7) + 
  labs(x = "Percent of Derived Alleles Decreasing Variability", 
       y = "Number of cis-eQTLs - Ctrl") + theme_classic()
Ctrl_cis_eQTL_FDI_plot

Ctrl_trans_eQTL_FDI_plot <- ggplot(Ctrl_trans_per_eQTL_FDI, aes(x = 100*(1-FDI))) + 
  annotate("text", y = Inf, x = 50, label = paste0("Total number of eQTLs = ",nrow(Ctrl_trans_per_eQTL_FDI)), size = 4, hjust = 0.5,vjust = 1)+
  geom_histogram(bins=5, fill = "lightblue", color = "black", alphsig = 0.7) + 
  labs(x = "Percent of Derived Alleles Decreasing Variability", 
       y = "Number of trans-eQTLs - Ctrl") + theme_classic()
Ctrl_trans_eQTL_FDI_plot

HS_cis_eQTL_FDI_plot <- ggplot(HS_cis_per_eQTL_FDI, aes(x = 100*(1-FDI))) + 
  annotate("text", y = Inf, x = 100, label = paste0("Total number of eQTLs = ",nrow(HS_cis_per_eQTL_FDI)), size = 4, hjust = 0.5,vjust = 1)+
  geom_histogram(bins=5, fill = "lightblue", color = "black", alphsig = 0.7) + 
  labs(x = "Percent of Derived Alleles Decreasing Variability", 
       y = "Number of cis-eQTLs - HS") + theme_classic()
HS_cis_eQTL_FDI_plot

HS_trans_eQTL_FDI_plot <- ggplot(HS_trans_per_eQTL_FDI, aes(x = 100*(1-FDI))) + 
  annotate("text", y = Inf, x = 50, label = paste0("Total number of eQTLs = ",nrow(HS_trans_per_eQTL_FDI)), size = 4, hjust = 0.5,vjust = 1)+
  geom_histogram(bins=5, fill = "lightblue", color = "black", alphsig = 0.7) + 
  labs(x = "Percent of Derived Alleles Decreasing Variability", 
       y = "Number of trans-eQTLs - HS") + theme_classic()
HS_trans_eQTL_FDI_plot
