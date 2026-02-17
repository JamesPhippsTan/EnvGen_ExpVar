# Concatenating the various downsamplings
# Running them together takes way too long
# Run them in parallel then merge

# Last Updated: 12/11/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(tidyr)
library(tibble)
library(matrixStats)
library(ggplot2)
library(ggrepel)
library(gamlss)

# Load gProfiler functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('gProfiler_Functions.R')


####################
##### Datasets #####
####################

# Background genes
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV")
background <- fread('AllGenes_DE_DV.csv',select = 1)
background_genes <- background$V1
number_of_expressed_genes = length(background_genes)

# Gene sets - take the final saved gene set at n=100
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_Downsampled")
load(file='DE_DV_Downsampled_Series_50100.RData')
gene_sets=list(canalised_only=subset(DVGsFinal,logFC.BCV<0 & padj.cpm >= 0.05),
               canalised_and_mean_change=subset(DVGsFinal,logFC.BCV<0 & padj.cpm <0.05))

gProfiler_results <- gProfiler_genesets(gene_sets = gene_sets, background_genes = background_genes)

# Note: canalised_only = 275; canalised_and_mean_change = 25'
# This has to be partially driven by lower numbers of DEGs

View(gene_sets[['canalised_only']])

#############################################
##### (1) GO Results - Full and Top 20 ######
#############################################

gProfiler_results <- gProfiler_genesets(gene_sets = gene_sets, background_genes = background_genes)
all_results <- gProfiler_results[[1]]
top_20_results <- gProfiler_results[[2]]

# Split dataframe into list of dataframes based on the dataset
all_results_split = split(all_results, f = all_results$dataset)
top_20_results_split = split(top_20_results, f = top_20_results$dataset)

#############################################
##### (2) Plotting Top 20 GO Results ########
#############################################

########################################
##### (2a) Top 20 GO:BP Results ########
########################################

GOBPplots <- lapply(names(top_20_results_split), 
                    plot_top_GO_refined, 
                    data = top_20_results_split, 
                    GO_type = "GO:BP")
names(GOBPplots) <- names(top_20_results_split)

GOBPplots[['canalised_only']] # neurogenesis
GOBPplots[['canalised_and_mean_change']] # nothing...

########################################
##### (2b) Top 20 GO:CC Results ########
########################################

GOCCplots <- lapply(names(top_20_results_split), 
                    plot_top_GO, 
                    data = top_20_results_split, 
                    GO_type = "GO:CC")
names(GOCCplots) <- names(top_20_results_split)

list(GOCCplots)


########################################
##### (2c) Top 20 GO:MF Results ########
########################################

GOMFplots <- lapply(names(top_20_results_split), 
                    plot_top_GO, 
                    data = top_20_results_split, 
                    GO_type = "GO:MF")
names(GOMFplots) <- names(top_20_results_split)

list(GOMFplots)


# Sure the metamorphosis enrichment go away but the neurogenesis ones stay
# What did I run the last time to erase the enrichment...
# Was the argument just that the N15 and outbred distributions of the genes were similar?