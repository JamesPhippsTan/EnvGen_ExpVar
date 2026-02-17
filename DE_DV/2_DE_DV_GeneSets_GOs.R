# Gene set enrichment analysis of various gene sets generated under the DE DV analysis

# Script by James Tan

# Last Updated: 2/12/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_GeneSets")
#load(file='DE_DV_GeneSets_GOs.RData')

library(gprofiler2)
library(ggplot2)
library(stringr)
library(rrvgo)
library(data.table)
library(svglite)
library(tidyverse)

# Load gProfiler functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('gProfiler_Functions.R')

######################
##### Datasets #######
######################

# Background genes
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV")
background <- fread('AllGenes_DE_DV.csv',select = 1)
background_genes <- background$V1
number_of_expressed_genes = length(background_genes)

# Gene sets
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_GeneSets")
gene_set_names <- list.files(pattern = "\\.csv$")
gene_sets <- lapply(gene_set_names, read.csv,row.names=1)
names(gene_sets) <- gsub("\\.csv$", "", gene_set_names)

# Include number of genes in the gene set names
gene_sets_summary <- do.call(rbind, lapply(names(gene_sets), function(item_name) {
  item_length <- nrow(gene_sets[[item_name]])  
  data.frame(Category = item_name, Number_of_genes = item_length)
}))
names(gene_sets) <- paste(gene_sets_summary$Category,"-",gene_sets_summary$Number_of_genes,"Genes")

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

setdiff(names(gene_sets),names(top_20_results_split))
# The large decanalized FC genes do not have any enrichments
# The large canalized FC genes do not have any enrichments

########################################
##### (2a) Top 20 GO:BP Results ########
########################################

GOBPplots <- lapply(names(top_20_results_split), 
                    plot_top_GO_refined, 
                    data = top_20_results_split, 
                    GO_type = "GO:BP")
names(GOBPplots) <- names(top_20_results_split)

# Which ones are interesting?
list(GOBPplots)


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

########################################
##### (2d) Top 20 KEGG Results #########
########################################

KEGGplots <- lapply(names(top_20_results_split), 
                    plot_top_GO, 
                    data = top_20_results_split, 
                    GO_type = "KEGG")
names(KEGGplots) <- names(top_20_results_split)

list(KEGGplots)

########################################
##### (2e) Top 20 TF Results ###########
########################################

TFplots <- lapply(names(top_20_results_split), 
                  plot_top_GO, 
                  data = top_20_results_split, 
                  GO_type = "TF")
names(TFplots) <- names(top_20_results_split)

list(TFplots)
# Increased - Trl motif
# Decreased (and decreased large) - mtTFA - mitochondrial transcription factor

########################################
##### (2f) Top 20 WP Results ###########
########################################

WPplots <- lapply(names(top_20_results_split), 
                  plot_top_GO, 
                  data = top_20_results_split, 
                  GO_type = "WP")
names(WPplots) <- names(top_20_results_split)

list(WPplots)

########################################
##### (2g) Top 20 miRNA Results ########
########################################

miRNAplots <- lapply(names(top_20_results_split), 
                     plot_top_GO, 
                     data = top_20_results_split, 
                     GO_type = "MIRNA")
names(miRNAplots) <- names(top_20_results_split)

list(miRNAplots)

View(all_results_split[["ADVGs_canalized - 789 Genes"]])
# How many genes from miRNA ROOT - 18 genes out of 70 genes

########################################
##### (2h) Top 20 HP Results ###########
########################################

HPplots <- lapply(names(top_20_results_split), 
                  plot_top_GO, 
                  data = top_20_results_split, 
                  GO_type = "HP")
names(HPplots) <- names(top_20_results_split)

list(HPplots)

#########################################
##### (3) 2 gene sets of interest #######
#########################################

# I want refined plots for the most interesting 2 gene sets
Canalized_genes_samemean_plot <- 
  plot_top_GO_refined(data=top_20_results_split, 
                      GO_type='GO:BP',
                      dataset_name="DVGs_canalized - 390 Genes", 
                      min_point_size = 1,
                      max_point_size = 5) 
Canalized_genes_samemean_plot

Canalized_genes_diffmean_plot <- 
  plot_top_GO_refined(data=top_20_results_split, 
                      GO_type='GO:BP',
                      dataset_name="DEDVGs_canalized - 399 Genes", 
                      min_point_size = 1,
                      max_point_size = 5) 
Canalized_genes_diffmean_plot


####################################
##### Save all the plots ###########
####################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_GeneSets")

# All gene sets
for (gene_set in names(top_20_results_split)){
  ggsave(plot=GOBPplots[[gene_set]],filename=paste0(gene_set,"_GOBPplots.png"),width=6,height = 4)
  ggsave(plot=GOCCplots[[gene_set]],filename=paste0(gene_set,"_GOCCplots.png"),width=6,height = 4)
  ggsave(plot=GOMFplots[[gene_set]],filename=paste0(gene_set,"_GOMFplots.png"),width=6,height = 4)
  ggsave(plot=TFplots[[gene_set]],filename=paste0(gene_set,"_TFplots.png"),width=6,height = 4)
  ggsave(plot=miRNAplots[[gene_set]],filename=paste0(gene_set,"_miRNAplots.png"),width=6,height = 4)
}

# Top 2 most interesting sets
ggsave(plot=Canalized_genes_samemean_plot,filename='Canalized_genes_samemean_plot.svg',width=5.5,height = 3,dpi=300) 
ggsave(plot=Canalized_genes_diffmean_plot,filename='Canalized_genes_diffmean_plot.svg',width=5.5,height = 3,dpi=300)

# Table of reported gene sets
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_GeneSets\\GO_Tables/")
DVGsets <- c("DEDVGs_canalized - 399 Genes","DVGs_canalized - 390 Genes","DVGs_decanalized - 3444 Genes")
for (DVG_set in DVGsets){
  sanitized_GO_table <- sanitize_for_csv(all_results_split[[DVG_set]])
  write.csv(sanitized_GO_table,file = paste0('GO_Table_',DVG_set,'.csv'),row.names = F)
}
save.image(file='DE_DV_GeneSets_GOs.RData')

############### Extras ###################

# Checking out the genes that increase in mean expression
DEGs_increased_table <- sanitize_for_csv(all_results_split[["DEGs_increased - 2766 Genes"]])
write.csv(DEGs_increased_table,file = paste0('GO_Table_DEGs_increased - 2766 Genes.csv'),row.names = F)

# Compress semantically similar terms for categories with >20 GO:BP terms
GOBP_sem_sim <- lapply(names(all_results_split), 
                       top_GO_sim_matrix, 
                       data = all_results_split, 
                       GO_label = "BP")
names(GOBP_sem_sim) <- names(all_results_split)
GOBP_sem_sim <-  GOBP_sem_sim[sapply(GOBP_sem_sim, length) >= 2]
GOBP_reduced <- lapply(names(GOBP_sem_sim), 
                       top_GO_reduced, 
                       simMatrix = GOBP_sem_sim,
                       similarity = 0.3,
                       GO_label = "BP") # get the reduced terms
names(GOBP_reduced) <- names(GOBP_sem_sim)
# Doesn't do too much, may as well just report top 20
