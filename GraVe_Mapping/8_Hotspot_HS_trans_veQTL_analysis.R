# Investigating top hotspot veQTL that form the tail of the HS trans-veQTL 

# Last Updated: 8/8/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls(all = T))

library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggridges)
library(reshape2)
library(scales)
library(svglite)


# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs\\HS_trans_veQTL_analyses/")
load(file = '8_Hotspot_HS_trans_veQTL_analysis.RData')

# Functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('gProfiler_Functions.R')
source('QTL_Analysis_Functions.R')
source('Variability_Functions.R')

#######################
##### Datasets ########
#######################

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
load(file = '5_Investigating_eQTL_and_veQTL.RData')

# Load hotspot SNPs
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs/")
hotspot_SNPs <- read.csv('hotspot_SNPs.csv')
hotspot_SNPs <- hotspot_SNPs[,-1]

# Significant eQTL and veQTL to get genes for each SNP
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
HS_trans_veQTL_sig <-read.csv("HS_trans_veQTL_sig.csv")
# Split the SNP names into SNP chromosomes and positions
HS_trans_veQTL_sig <-extract_chromosome_and_pos(HS_trans_veQTL_sig,input_col = 'SNP',keep_original = T)

# TFBSs and coordinates
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
All_genes <-fread("Gene_MeanVar_Table.csv",select = 1)

#######################################
#### (1) Select top HS trans-veQTL ####
#######################################

HS_trans_veQTL_df <- na.omit(data.frame(row.names = hotspot_SNPs$SNP,Number_of_genes=hotspot_SNPs$HS_trans_veQTLs))

summary(HS_trans_veQTL_df$Number_of_genes)

# Replot with the 4 cutoffs
binwidth <- ifelse(max(HS_trans_veQTL_df$Number_of_genes)<10,1,
                   (max(HS_trans_veQTL_df$Number_of_genes) - min(HS_trans_veQTL_df$Number_of_genes)) / max(HS_trans_veQTL_df$Number_of_genes))
veQTL_pleiotropy <- ggplot(HS_trans_veQTL_df) +
  geom_histogram(
    aes(x = Number_of_genes),
    bins = 10,
    binwidth = binwidth,
    center = 0,
    position = "identity",
    color = "black",
    fill= '#611BB8'
  ) +
  theme_classic() +
  labs(
    title = 'HS trans-veQTL',
    subtitle = paste0(nrow(HS_trans_veQTL_df),' SNPs'),
    x = "# of regulated genes per SNP",
    y = "# of SNPs (log scale)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )+
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10),
    breaks = c(0, 1, 10, 100, 1000, 10000)
  )+geom_vline(xintercept = c(10,40,100,200),color='red', linetype = 'dashed')

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs")
ggsave(plot = veQTL_pleiotropy,filename = paste0('trans-veQTL_MAF_all_labelled.svg'),dpi=300,width = 4,height = 2)

# Let us run a gradient of 4 sizes
HS_trans_veQTL_10 <- subset(HS_trans_veQTL_df,Number_of_genes>10)
# 10K SNPs more than 10 genes...
# Let us try 40 out, given this is the maximum number of eQTL
HS_trans_veQTL_40 <- subset(HS_trans_veQTL_df,Number_of_genes>40)
# 1735 SNPs more than 40 genes...
HS_trans_veQTL_100 <- subset(HS_trans_veQTL_df,Number_of_genes>100)
# 300 SNPs more than 100 genes...
HS_trans_veQTL_200 <- subset(HS_trans_veQTL_df,Number_of_genes>200)
# 50 SNPs more than 100 genes...

# Maybe they tag some weird region?
View(data.frame(X=rownames(HS_trans_veQTL_100)))

# Make them into a list
veQTL_sets <- list(HS_trans_veQTL_10=rownames(HS_trans_veQTL_10),
            HS_trans_veQTL_40=rownames(HS_trans_veQTL_40),
            HS_trans_veQTL_100=rownames(HS_trans_veQTL_100),
            HS_trans_veQTL_200=rownames(HS_trans_veQTL_200))

# Let us look at their MAF distributions
for (i in seq_along(veQTL_sets)) {
  set_name <- names(veQTL_sets)[i]
  set <- veQTL_sets[[i]]
  title <- paste0('>',gsub('HS_trans_veQTL_','',set_name),' genes','\n',
                  length(set),' SNPs')
  
  max_count <- max(hist(MAF[set, ]$HS_MAF, breaks = 10, plot = FALSE)$counts)
  print(ks.test(MAF[set, ]$HS_MAF,MAF[rownames(HS_trans_veQTL_df), ]$HS_MAF))
  MAF_plot <- ggplot(data = MAF[set, ], aes(x = HS_MAF)) +
    geom_histogram(bins = 10, fill = "#611BB8", color = "black") +
    labs(x = 'MAF', y = 'Count',size=15) +
    annotate("text",
             x = 0.5,                           # Adjust if needed based on range of HS_MAF
             y = max_count,                   # Use estimated max count for top y
             label = title,
             size = 4,
             hjust = 1, vjust = 1) +          # Top-right alignment
    theme_classic()
  
  setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs\\HS_trans_veQTL_analyses")
  ggsave(plot = MAF_plot,filename = paste0('trans-veQTL_MAF_',set_name,'.svg'),dpi=300,width = 2,height = 2)
}
# They are not biased towards the lower end...even when considering the top 50
# They are rather evenly distributed across the MAF sizes
# It is not purely an artefact of imbalanced sample size genotypes being more susceptible to being called a veQTL

##################################################
#### (2) Gene set enrichment for top 50 veQTL ####
##################################################

# Let us then focus on this top set of 50 for now...
# We can run the following analysis to get 50 gene set enrichments...
top_hotspot_SNPs <- rownames(HS_trans_veQTL_200)

# Make a list of gene sets for each veQTL - there should be 50 in total
gene_sets <- list()
for (hotspot_SNP in top_hotspot_SNPs){
  gene_sets[[hotspot_SNP]] <- subset(HS_trans_veQTL_sig,SNP == hotspot_SNP)$GENE 
}

# Get summary of gene sets
gene_sets_summary <- do.call(rbind, lapply(names(gene_sets), function(item_name) {
  item_length <- length(gene_sets[[item_name]])  # Get the length of the item
  data.frame(Category = item_name, Number_of_genes = item_length)  # Create a data frame
}))

# Add the number of genes to the gene set name
names(gene_sets) <- paste(gene_sets_summary$Category,"-",gene_sets_summary$Number_of_genes," vGenes")

# Get the background genes - all expressed genes in the head
background_genes <- All_genes$phenotype_id
number_of_expressed_genes = length(background_genes)

# Run gProfiler
gProfiler_results <- gProfiler_genesets(gene_sets = gene_sets, background_genes = background_genes)
all_results <- gProfiler_results[[1]]
top_20_results <- gProfiler_results[[2]]

# Split dataframe into list of dataframes based on the dataset
all_results_split = split(all_results, f = all_results$dataset)
top_20_results_split = split(top_20_results, f = top_20_results$dataset)
unique(all_results$source)

# GO:BPs
GOBPplots <- lapply(names(top_20_results_split), 
                    plot_top_GO, 
                    data = top_20_results_split, 
                    GO_type = "GO:BP")
names(GOBPplots) <- names(top_20_results_split)
list(GOBPplots) 

# Save enrichments
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs\\HS_trans_veQTL_analyses/")
write.csv(gene_sets_summary,'Hotspot_SNPs_gene_sets_summary.csv')
for (gene_set in names(top_20_results_split)){
  ggsave(plot=GOBPplots[[gene_set]],filename=paste0("GOBPplots_",as.character(gene_set),".svg"),dpi=300,,width=6,height = 4)
}
# BPs- mostly RNA processing and transcription
# How can this be? 50 genes, 40 of them tell me the gene sets are transcription
# TFs - SRYDELTA, DREF, BEAF32, pnr
# Serendipity and oogenesis - bicoid controller 
# DNA replication factor - DNAP
# Boundary element-associated factor of 32kD (insulator protein!) - terminal neuron dev, hippo, eye development
# https://flybase.org/reports/FBrf0232800 
# Pannier (proneural TF) - terminal neuron dev
# Maybe the genes are similar...
# Just the BPs
all_results_BP <- subset(all_results,source=='GO:BP')
View(all_results_BP)
write.csv(x=sanitize_for_csv(all_results_BP),file = 'pleiotrophicveQTL10_all_results_BP.csv')

# Let us quantify the most common GO:BP terms
BP_term_counts <- as.data.frame(table(all_results_BP$term_name))
colnames(BP_term_counts) <- c("BP_term_name", "count")
View(BP_term_counts)

# Terms present in 10 or more SNPs
BP_term_counts_10 <- subset(BP_term_counts,count>10)
BP_term_counts_10 <- BP_term_counts_10[order(BP_term_counts_10$count),]
BP_term_counts_10$percent=100*(BP_term_counts_10$count/55)
BP_term_counts_10$BP_term_name <- factor(
  BP_term_counts_10$BP_term_name,
  levels = BP_term_counts_10$BP_term_name[order(BP_term_counts_10$count,decreasing = T)]
)
BP_10more_plot <- ggplot(BP_term_counts_10, aes(x = BP_term_name, y = percent)) +
  geom_col(fill = "#611BB8") +  # Use geom_col() for pre-computed counts
  labs(x = "15 most common GO:BP terms", y = "Percent out of 55 SNPs") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)  # Angled at 45°
  )
ggsave(plot=BP_10more_plot,filename='BP10ormoregenesets.svg',dpi=300,,width=7,height = 7)

####################################################################
#### (3) Are the genes regulated by pleiotropic veQTL the same? ####
####################################################################

# Worth doing an all against all overlap of the gene sets to determine if these are the same sets of genes
# That tend to have veQTL in HS conditions

# Function to compute Jaccard index between two sets
jaccard_index <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
}

# Initialize an empty matrix
n <- length(gene_sets)
similarity_matrix <- matrix(NA, nrow = n, ncol = n)

# Fill in the matrix
for (i in 1:n) {
  for (j in 1:n) {
    similarity_matrix[i, j] <- jaccard_index(gene_sets[[i]], gene_sets[[j]])
  }
}

# Add row and column names for clarity
names_list <- names(gene_sets)
rownames(similarity_matrix) <- names_list
colnames(similarity_matrix) <- names_list

# View the matrix
View(similarity_matrix)
max_non_one <- max(similarity_matrix[similarity_matrix != 1], na.rm = TRUE)
max_non_one

melted_sim <- melt(similarity_matrix)
colnames(melted_sim) <- c("List1", "List2", "Similarity")

melted_sim <- melted_sim %>%
  mutate(
    List1 = factor(List1, levels = unique(c(List1, List2))),
    List2 = factor(List2, levels = unique(c(List1, List2))))
    
# Keep upper triangle (List1 <= List2 to include diagonal)
melted_sim_upper <- melted_sim %>%
      filter(as.numeric(List1) <= as.numeric(List2))
    
# Plot
Jaccard_similarity_plot <- ggplot(
      melted_sim_upper,
      aes(x = List1, y = List2, fill = Similarity)
    ) +
      geom_tile(color = "grey80") +
      scale_fill_gradient(
        low = "white", 
        high = "blue", 
        limits = c(0, 1),
        na.value = "white"  # Hide NA values (unnecessary here but safe)
      ) +
      geom_text(
        aes(label = round(Similarity, 2)), 
        size = 1.5, 
        color = "black"  # Ensure text is visible on light tiles
      ) +
      theme_minimal(base_size = 6) +
      theme(
        axis.text.x = element_text(
          angle = 90, 
          vjust = 0.5, 
          hjust = 1, 
          size = 4
        ),
        axis.text.y = element_text(size = 4),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right"
      ) +
      labs(
        fill = "Similarity"
      ) +
      coord_fixed() + scale_x_discrete(position = "top")
    
# Print plot
print(Jaccard_similarity_plot)
ggsave(plot=Jaccard_similarity_plot,filename='Jaccard_similarity_plot.svg',dpi=300,width=10,height = 7.5)

####################################################################
#### (4) Which genes are the most recurrent within this set? #######
####################################################################

all_genes <- unlist(gene_sets)

# For unique counts per gene (across sets, not total mentions):
# Make a logical matrix: rows = genes, columns = sets
all_unique_genes <- unique(all_genes)
gene_in_set_matrix <- sapply(gene_sets, function(set) all_unique_genes %in% set)

# Now count how many sets each gene appears in
gene_set_counts <- rowSums(gene_in_set_matrix)

# Convert to a named vector (or data frame)
names(gene_set_counts) <- all_unique_genes
gene_set_counts_df <- data.frame(Gene = all_unique_genes, Count = gene_set_counts)

# View the genes that appear the most often
View(gene_set_counts_df)
# The highest genes appear in 40 sets!
# These are genes that have many HS vSNPs
# This corresponds to the third quartile of the number of SNPs per gene
# (A few thousand is the maximum)

# Flybase them
# Taspase1 - trithorax endopeptidase expression regulator - G- bacteria
# Fibrinogen

##############################################
##### Save the working environment ###########
##############################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs\\HS_trans_veQTL_analyses/")
save.image(file='8_Hotspot_HS_trans_veQTL_analysis.RData')


# Look at genes regulated by ecdysone-kinase-22 hotspot region
ecdy22kin <- hotspot_SNPs[grepl("3R252", hotspot_SNPs$SNP) & hotspot_SNPs$Ctrl_cis_eQTLs > 6, ]
View(ecdy22kin)

ecdy22kin_SNPs <- as.vector(na.omit(ecdy22kin$SNP))
View(Ctrl_cis_eqtl_sig)

# Control cis
Ctrl_cisGenes <- lapply(ecdy22kin_SNPs, function(s) Ctrl_cis_eqtl_sig$phenotype_id[Ctrl_cis_eqtl_sig$SNP == s])
names(Ctrl_cisGenes) <- ecdy22kin_SNPs
View(Ctrl_cisGenes)
Ctrl_cisGenes_unique <- unique(unlist(Ctrl_cisGenes))

# HS cis
HS_cisGenes <- lapply(ecdy22kin_SNPs, function(s) HS_cis_eqtl_sig$phenotype_id[HS_cis_eqtl_sig$SNP == s])
names(HS_cisGenes) <- ecdy22kin_SNPs
View(HS_cisGenes)
HS_cisGenes_unique <- unique(unlist(HS_cisGenes))

# Overlap
Ctrl_cisGenes_unique # 12 genes
HS_cisGenes_unique # 11 genes
intersect(HS_cisGenes_unique,Ctrl_cisGenes_unique) # 10 genes
# Dro2-0, Dro5-1, Dro5-6, Dro4-0, Dro5-5
# Dro11-0, Dro13-1, Dro5-8, Dro5-7, Dro5-2
setdiff(HS_cisGenes_unique,Ctrl_cisGenes_unique) # Dro10-0
setdiff(Ctrl_cisGenes_unique,HS_cisGenes_unique) # Dro3-0 and Dro9-0
# Note: none of these are developmentally lethal in RNAi lines
# One affects the number of closing adults