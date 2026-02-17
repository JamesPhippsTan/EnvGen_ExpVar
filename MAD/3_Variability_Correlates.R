# Transcriptome variance analyses of Adult Netherlands D.mel heads on normal vs high sugar diets
# Identifying correlations between various metrics of variability and gene-specific features

# Script by James Tan

# Last Updated: 15/9/25

################################
##### Packages and Setup #######
################################

rm(list = ls())

library(tidyr)
library(dplyr)
library(phylomapr)
library(ggplot2)
library(ppcor)
library(tibble)
library(STRINGdb)
library(igraph)
library(orthologr)
library(biomartr)
library(rdiamond)
library(DescTools)
library(grid)
library(gridExtra)
library(EnvStats)
library(patchwork)
library(ggsignif)
library(rlang)
library(svglite)


# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD\\")
load(file='Variability_Correlates.RData')

# Load functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source("Quantile_Functions.R")
source("Variability_Functions.R")

######################
##### Datasets #######
######################

# Use 'gene_ID' for merging

# Per-gene means and variabilities
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
Gene_MeanVar_Table <- read.csv("Gene_MeanVar_Table.csv",row.names = 1)
Gene_MeanVar_Table <- rownames_to_column(Gene_MeanVar_Table,var='gene_ID')

# DE DV test
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV")
DE_DV_Table <- read.csv("AllGenes_DE_DV.csv",row.names = 1)
DE_DV_Table <- rownames_to_column(DE_DV_Table,var='gene_ID')

# Get variability correlates
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD\\Variability Correlates Data")

# Map of flybase IDs
Fbpp_Fbtr_Fbgn_Map <- read.csv('flybase_fbgn_fbtr_fbpp_expanded_fb_2024_03.csv') 

# (1) Gene length, heritability, nucleotide diversity - 'gene'
eQTL_covars <- read.csv("eQTL_covariates_head.csv")
names(eQTL_covars)[names(eQTL_covars) == "gene"] <- "gene_ID"

# (2) Number of transcripts - 'gene_ID'
Transcript_number_per_gene <- function(gene_df){
  split_gene_df <- split(gene_df,f=gene_df[,'gene_ID'])
  transcripts_per_gene <- as.data.frame(do.call(rbind,lapply(split_gene_df, nrow)))
  transcripts_per_gene$gene_ID <- rownames(transcripts_per_gene)
  colnames(transcripts_per_gene)[1] <- 'number_of_unique_transcripts'
  return(transcripts_per_gene)
}
Transcript_number <- Transcript_number_per_gene(Fbpp_Fbtr_Fbgn_Map)

# (3) Essentiality - 'locus'
Essentiality <- read.table("gene_essentiality.txt",header=T,fill = TRUE,sep='\t')
Essentiality <- subset(Essentiality, taxaID == 7227)
Essentiality$essential <- ifelse(Essentiality$essentiality=='E',1,0)
Essentiality$conditional_essential <- ifelse(Essentiality$essentiality=='C',1,0)
Essentiality$non_essential <- ifelse(Essentiality$essentiality=='NE',1,0)
names(Essentiality)[names(Essentiality) == "locus"] <- "gene_ID"

# (4) dNdS - 'transcript_id'
# Run diamond
#Dmel_filepath <- biomartr::getCDS( db = "ensembl", organism = "Drosophila melanogaster")
#Dpse_filepath <- biomartr::getCDS( db = "ensembl", organism = "Drosophila pseudoobscura")
#setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD\\Variance Correlates Data")
#dNdS <-
#  dNdS(
#    query_file      = Dmel_filepath,
#    subject_file    = Dpse_filepath,
#    delete_corrupt_cds = TRUE,
#    ortho_detection = "RBH", # perform DIAMOND best reciprocal hit orthology inference
#    aa_aln_type     = "pairwise", # perform pairwise global alignments of AA seqs
#    aa_aln_tool     = "NW", # using Needleman-Wunsch
#    codon_aln_tool  = "pal2nal", # perform codon alignments using the tool Pal2Nal
#    dnds_est.method = "Comeron", # use Comeron's method for dN/dS inference
#    store_locally = T # store the alignment locally) 
# Note - I copied the diamond.exe file into my working directory
dNdS <- read.csv('dNdS.csv') 
colnames(dNdS)[1] <- 'transcript_ID'
dNdS <- merge(dNdS[,c('transcript_ID','dNdS')],Fbpp_Fbtr_Fbgn_Map[,c('transcript_ID','gene_ID')],by='transcript_ID')
dNdS <- distinct(dNdS, gene_ID, .keep_all = TRUE)

# (4) GenEra-inferred gene ages - 'gene_ID'
# (mergings suggested by Josue to have a balance between sample size and taxonomic distinction)
Dmel_PhyloMap_GenEra_ENSEMBL <- phylomapr::Drosophila_melanogaster_ENSEMBL_BDGP5.PhyloMap
colnames(Dmel_PhyloMap_GenEra_ENSEMBL)[2] <- 'polypeptide_ID'
Dmel_PhyloMap_GenEra <- merge(Dmel_PhyloMap_GenEra_ENSEMBL[,c('polypeptide_ID','Phylostratum')],Fbpp_Fbtr_Fbgn_Map[,c('polypeptide_ID','gene_ID')],by='polypeptide_ID')
# 1 - LUCA; 2 - Eukarya; 3 - Opisthokonts; 4 - Metazoa
# 5 - Bilateria; 6 - Protostomia; 7 - Ecdysozoa
Dmel_PhyloMap_GenEra$Phylostratum <- replace(Dmel_PhyloMap_GenEra$Phylostratum, Dmel_PhyloMap_GenEra$Phylostratum > 7 & Dmel_PhyloMap_GenEra$Phylostratum < 12, 8)
# 8 - Panarthropoda - Hexapoda (i.e., Insecta)
Dmel_PhyloMap_GenEra$Phylostratum <- replace(Dmel_PhyloMap_GenEra$Phylostratum, Dmel_PhyloMap_GenEra$Phylostratum == 12, 9)
# 9 - Pterygota
Dmel_PhyloMap_GenEra$Phylostratum <- replace(Dmel_PhyloMap_GenEra$Phylostratum, Dmel_PhyloMap_GenEra$Phylostratum > 12 & Dmel_PhyloMap_GenEra$Phylostratum < 16, 10)
# 10 - Neoptera and Endopterygota
Dmel_PhyloMap_GenEra$Phylostratum <- replace(Dmel_PhyloMap_GenEra$Phylostratum, Dmel_PhyloMap_GenEra$Phylostratum == 16, 11)
# 11 - Diptera
Dmel_PhyloMap_GenEra$Phylostratum <- replace(Dmel_PhyloMap_GenEra$Phylostratum, Dmel_PhyloMap_GenEra$Phylostratum > 16 & Dmel_PhyloMap_GenEra$Phylostratum < 22, 12)
# 12 - Brachycera - Acalyptratae
Dmel_PhyloMap_GenEra$Phylostratum <- replace(Dmel_PhyloMap_GenEra$Phylostratum, Dmel_PhyloMap_GenEra$Phylostratum == 22, 13)
# 13 Drosophilinae
Dmel_PhyloMap_GenEra$Phylostratum <- replace(Dmel_PhyloMap_GenEra$Phylostratum, Dmel_PhyloMap_GenEra$Phylostratum > 22 & Dmel_PhyloMap_GenEra$Phylostratum < 27, 14)
# 14 Drosophila genus - subgroup
Dmel_PhyloMap_GenEra$Phylostratum <- replace(Dmel_PhyloMap_GenEra$Phylostratum, Dmel_PhyloMap_GenEra$Phylostratum >= 27, 15)
# 15 Drosophila melanogaster (species)
Dmel_PhyloMap_GenEra$Phylostratum <- as.factor(Dmel_PhyloMap_GenEra$Phylostratum)
Dmel_PhyloMap_GenEra <- distinct(Dmel_PhyloMap_GenEra, gene_ID, .keep_all = TRUE)



# (5) Protein abundance - 'gene_ID'
abundances <- read.csv('Proteome_Main.csv',sep=',') 
mapped_IDs_1 <- read.table('proteome_1_2025_01_21.tsv',header = T) 
mapped_IDs_2 <- read.table('proteome_2_2025_01_21.tsv',header = T) 
mapped_abundances_1 <- merge(abundances,mapped_IDs_1,by='UniProt_ID',all=T)
mapped_abundances_2 <- merge(mapped_abundances_1,mapped_IDs_2,by='UniProt_Acc',all=T)
mapped_abundances_3 <- mapped_abundances_2 %>%
  mutate(gene_ID = coalesce(Gene_name_1, Gene_name_2)) 
mapped_abundances_final <- mapped_abundances_3 %>%
  group_by(gene_ID) %>%
  summarise(max_protein_level = max(emPAI),min_protein_level = min(emPAI),protein_level_range = max(emPAI)-min(emPAI), .groups = "drop")

# (6) Number of PPIs from StringDB - 'locus'
string_db <- STRINGdb$new(species=7227,version='12',score_threshold=400,network_type='physical')
dros_graph <- string_db$get_graph()
V(dros_graph)$degree <- degree(dros_graph) # Degree - how many PPIs
V(dros_graph)$hubs <- hub_score(dros_graph)$vector               
V(dros_graph)$closeness <- closeness(dros_graph)                
V(dros_graph)$betweenness <- betweenness(dros_graph)          
V(dros_graph)$transitivity <- transitivity(dros_graph,type='local')
PPI <- data.frame(polypeptide_ID   = gsub('7227.','',V(dros_graph)$name),
                         PPI_degree      = V(dros_graph)$degree,
                         PPI_closeness   = V(dros_graph)$closeness,
                         PPI_betweenness = V(dros_graph)$betweenness,
                         PPI_hubness = V(dros_graph)$hubs,
                         PPI_transitivity = V(dros_graph)$transitivity)
PPI <- merge(PPI,Fbpp_Fbtr_Fbgn_Map[,c('polypeptide_ID','gene_ID')],by='polypeptide_ID')

# (7) Molecular network connectivity - "FBgn"
network_connectivity <- read.table("pleio_putatively_adaptive_gene.txt",sep = "\t",header = T)
names(network_connectivity)[names(network_connectivity) == "FBgn"] <- "gene_ID"
network_connectivity$transcriptional_connectivity <- network_connectivity$connectivity

# (8) Genomic features (TATA Box, Promoter Shape, etc.) - 'gene_id'
Genomic_features_full <- read.table("master_table.gene_radius_dhs.prox500.dist10000.txt",header = T)
# We want only features that are gene-specific and independent of stage
# Identify columns where all values are identical within each group
check_identical <- function(col, group) {tapply(col, group, function(x) length(unique(x)) == 1)}
identical_columns <- sapply(Genomic_features_full[,-1], function(col) all(check_identical(col, Genomic_features_full$gene_id)))
# Subset dataframe to keep only these columns
subset_df <- Genomic_features_full[, c(TRUE, identical_columns)]
# Eliminate repeat rows using aggregate, which takes the first row of each gene
Genomic_features <- aggregate(. ~ gene_id, data = subset_df, FUN = function(x) x[1])
names(Genomic_features)[names(Genomic_features) == "gene_id"] <- "gene_ID"

#####################################
##### (1) Create mega-dataset #######
#####################################

# List of all dataframes to merge
Variability_correlates_list <- list(Gene_MeanVar_Table, 
                   eQTL_covars,
                   Transcript_number,
                   Essentiality,
                   DE_DV_Table,
                   dNdS,
                   Dmel_PhyloMap_GenEra,
                   PPI, 
                   tau, 
                   network_connectivity,
                   mapped_abundances_final,
                   Genomic_features)
Variability_correlates_df <- Reduce(function(x, y) merge(x, y, by = "gene_ID", all.x = TRUE), Variability_correlates_list)

#######################################################
##### (2) Select mean and variability metric ##########
#######################################################

# Decide on mean and variability metric
var_metric <- "VST_MAD"
mean_metric <- "VST_Mean"

# For plot saving
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD\\")

###########################################################################
##### (2) Numeric variables - Linear plots and rank correlations ##########
###########################################################################

Numeric_correlates <- c("gene_length",
                        "averagepi",
                        "number_of_unique_transcripts",
                        "major_shape_ind",
                     "tss_gc",
                     "gene_gc",
                     "num_dhs_conditions.prox",
                     "dNdS",
                     "PPI_degree",
                     "PPI_hubness",
                     "PPI_closeness",
                     "PPI_betweenness",
                     "PPI_transitivity",
                     "transcriptional_connectivity",
                     "h2_h")

X_axis_names <- c("Gene length / Rank",
                        "Nuceleotide diversity (pi) / Rank",
                        "Number of unique transcripts / Rank",
                        "Major promoter shape index / Rank",
                        "TSS GC content / Rank",
                        "Genic GC content / Rank",
                        "Number of conditions with detected DHS / Rank",
                        "dN/dS / Rank",
                        "PPI Network Degree / Rank",
                        "PPI Network Number of Hubs / Rank",
                        "PPI Network Closeness / Rank",
                        "PPI Network Betweenness / Rank",
                        "PPI Network Transitivity / Rank",
                        "Transcriptional Network Connectivity / Rank",
                  "Heritability / Rank")

Numeric_correlation_table <- data.frame(
  Condition = character(),
  Correlate = character(),
  Number_of_genes = numeric(),
  Spearman_Coefficient = numeric(),
  Spearman_P_Value = numeric(),
  Var_Spearman_Coefficient = numeric(),
  Var_Spearman_P_Value = numeric(),
  Mean_Spearman_Coefficient = numeric(),
  Mean_Spearman_P_value = numeric(),
  stringsAsFactors = FALSE
)

for (correlate_index in 1:length(Numeric_correlates)){
  
  correlate <- Numeric_correlates[correlate_index]
  meanplot_list <- list()
  varplot_list <- list()
  plot_list <- list()
  for (condition in c("Ctrl","HS")){
    
  var_category <- paste0(condition,"_",var_metric)
  mean_category <- paste0(condition,"_",mean_metric)
  
  df_subset <- Variability_correlates_df %>%
    dplyr::select(!!sym(var_category),!!sym(mean_category), !!sym(correlate)) %>%
    na.omit()
  
  df_subset[[var_category]] <- as.numeric(df_subset[[var_category]])
  df_subset[[mean_category]] <- as.numeric(df_subset[[mean_category]])
  df_subset[[correlate]] <- as.numeric(df_subset[[correlate]])
  correlate_quantile <- paste0(correlate," Quantile")
  df_subset <- make_quantiles(data=df_subset,
                              quantile_column = df_subset[[correlate]],
                                                  n_quantiles = 10,
                                                  name_of_quantiles = correlate_quantile)
  
  # Extract Spearman correlation coefficient and p-value
  spearman_test <- cor.test(y=df_subset[[var_category]], x=df_subset[[correlate]], method = "spearman")
  spearman_coef <- spearman_test$estimate
  spearman_p_value <- spearman_test$p.value

  # Extract partial correlation coefficient and p-value
  partial_test <- pcor.test(y=df_subset[[var_category]], x=df_subset[[correlate]], z=df_subset[[mean_category]], method = "spearman")
  partial_coef <- partial_test$estimate
  partial_p_value <- partial_test$p.value
  
  # Repeat spearman for mean
  mean_spearman_test <- cor.test(y=df_subset[[mean_category]], x=df_subset[[correlate]], method = "spearman")
  mean_spearman_coef <- mean_spearman_test$estimate
  mean_spearman_p_value <- mean_spearman_test$p.value
  
  # Add the results to the correlation table
  Numeric_correlation_table <- rbind(Numeric_correlation_table, data.frame(
    Condition = condition,
    Correlate = correlate,
    Number_of_genes = nrow(df_subset),
    Spearman_Coefficient = spearman_coef,
    Spearman_P_Value = spearman_p_value,
    Var_Spearman_Coefficient = partial_coef,
    Var_Spearman_P_Value = partial_p_value,
    Mean_Spearman_Coefficient = mean_spearman_coef,
    Mean_Spearman_P_value = mean_spearman_p_value,
    stringsAsFactors = FALSE
  ))
  
  # Define text
  Y_axis_name <- paste0(condition,' transcript level variability (MAD)')
  Y_axis_name_mean <- paste0(condition,' mean transcript level')
  X_axis_name <- X_axis_names[correlate_index]
  number_of_genes <- nrow(df_subset)
  stat_test_result <- paste0("Partial Spearman rho = ",format_statistic(partial_coef),"\n",
         "p = ",format_statistic(partial_p_value),"\n",'Number of genes = ',number_of_genes)
  mean_stat_test_result <- paste0("Spearman rho = ",format_statistic(mean_spearman_coef),"\n",
                             "p = ",format_statistic(mean_spearman_p_value),"\n",'Number of genes = ',number_of_genes)
  
  # Plot
  var_plot <- ggplot(df_subset, aes(y = !!sym(var_category), x = rank(!!sym(correlate)))) +
    geom_point(color='darkgrey',alpha=0.8) +
    geom_smooth(method = "lm", col = "blue",linetype='dashed')+
    ylab(Y_axis_name)+
    xlab(X_axis_name)+
    annotate("text", y = max(df_subset[[var_category]]), x = max(rank(df_subset[[correlate]])), label = stat_test_result, size = 2.35, hjust = 1,vjust = 1)+
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 5))+
    theme_classic()  
  
  mean_plot <- ggplot(df_subset, aes(y = !!sym(mean_category), x = rank(!!sym(correlate)))) +
    geom_point(color='darkgrey',alpha=0.8) +
    geom_smooth(method = "lm", col = "blue",linetype='dashed')+
    ylab(Y_axis_name_mean)+
    xlab(X_axis_name)+
    annotate("text", y = max(df_subset[[mean_category]]), x = max(rank(df_subset[[correlate]])), label = mean_stat_test_result, size = 2.35, hjust = 1,vjust = 1)+
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 5))+
    theme_classic()  

  plot_list[[mean_category]] <- mean_plot
  plot_list[[var_category]] <- var_plot
  meanplot_list[[mean_category]] <- mean_plot
  varplot_list[[var_category]] <- var_plot
  }
  # Save the plots
  mean_plots <- wrap_plots(meanplot_list, nrow = 2)  # Arrange in 2 columns
  print(mean_plots)
  ggsave(plot = mean_plots,filename = paste0("MeanCor_",correlate,"plot.svg"),width=3,height=6,dpi=300)
  var_plots <- wrap_plots(varplot_list, nrow = 2)  # Arrange in 2 columns
  print(var_plots)
  ggsave(plot = var_plots,filename = paste0("VarCor_",correlate,"plot.svg"),width=3,height=6,dpi=300)
  combined_plots <- wrap_plots(plot_list, nrow = 2)  # Arrange in 2 columns
}


#################################################################
##### (3) Categorical variables - Boxplots and Wilcox tests #####
#################################################################

# Table of values and medians
Categorical_test_table <- data.frame(
  Condition = character(),
  Variable = character(),
  Number_of_genes = numeric(),
  test_type = character(),
  Statistic = numeric(),
  P_Value = numeric(),
  Mean_Statistic = numeric(),
  Mean_P_Value = numeric(),
  Median_Group0 = numeric(),
  Median_Group1 = numeric(),
  Mean_Median_Group0 = numeric(),
  Mean_Median_Group1 = numeric(),
  stringsAsFactors = FALSE
)

# Split between the 2-factor categories and phylostratum 

# Binary variables
binary_vars <- c("ohler_maj.TATA", "is_housekeeping", "is_ubiquitous", "is_maternal")
X_axis_names <- c("Has TATA Box", "Housekeeping", "Ubiquitously expressed", "Maternally expressed")
format_statistic <- function(x, digits = 2) {
  if (is.null(x)) return("NA")
  if (x < 0.001) return(formatC(x, format = "e", digits = digits))
  return(round(x, digits))
}

color_list <- c("#E2CBD3", "#C7CCE8", "#C7CCE8", "#C7CCE8")

# Loop over binary variables
for (variable_index in seq_along(binary_vars)) {
  variable <- binary_vars[[variable_index]]
  meanplot_list <- list()
  varplot_list <- list()
  
  for (condition in c("Ctrl", "HS")) {
    var_category <- paste0(condition, "_", var_metric)
    mean_category <- paste0(condition, "_", mean_metric)
    
    df_subset <- Variability_correlates_df %>%
      dplyr::select(!!sym(var_category), !!sym(mean_category), !!sym(variable)) %>%
      na.omit()
    
    df_subset[[variable]] <- factor(ifelse(df_subset[[variable]] == 1, "Yes", "No"), levels = c("No", "Yes"))
    df_subset[[var_category]] <- as.numeric(df_subset[[var_category]])
    df_subset[[mean_category]] <- as.numeric(df_subset[[mean_category]])
    
    # Wilcoxon tests
    test <- wilcox.test(df_subset[[var_category]] ~ df_subset[[variable]])
    mean_test <- wilcox.test(df_subset[[mean_category]] ~ df_subset[[variable]])
    
    # Median values
    medians <- tapply(df_subset[[var_category]], df_subset[[variable]], median)
    mean_medians <- tapply(df_subset[[mean_category]], df_subset[[variable]], median)
    delta_median <- diff(medians)
    delta_mean_median <- diff(mean_medians)
    
    # Append to results table
    Categorical_test_table <- rbind(Categorical_test_table, data.frame(
      Condition = condition,
      Variable = variable,
      Number_of_genes = nrow(df_subset),
      test_type = "Wilcoxon W",
      Statistic = test$statistic,
      P_Value = test$p.value,
      Mean_Statistic = mean_test$statistic,
      Mean_P_Value = mean_test$p.value,
      Median_Group0 = medians[["No"]],
      Median_Group1 = medians[["Yes"]],
      Mean_Median_Group0 = mean_medians[["No"]],
      Mean_Median_Group1 = mean_medians[["Yes"]],
      stringsAsFactors = FALSE
    ))
    
    # Labels
    X_axis_name <- X_axis_names[[variable_index]]
    sig_label <- paste0("Δmedian = ", format_statistic(delta_median), "\np = ", format_statistic(test$p.value))
    mean_sig_label <- paste0("Δmedian = ", format_statistic(delta_mean_median), "\np = ", format_statistic(mean_test$p.value))
    
    # Group sizes
    group_counts <- table(df_subset[[variable]])
    group_names <- names(group_counts)
    count_labels <- c("No", "Yes")  # Clean labels
    
    # Fixed y-axis position for n labels (just below the boxplot area)
    y_text_var <- min(log10(df_subset[[var_category]]), na.rm = TRUE) - 0.1
    y_text_mean <- min(df_subset[[mean_category]], na.rm = TRUE) - 0.1
    
    # Position for significance annotation (above plot)
    y_position_var <- max(log10(df_subset[[var_category]]), na.rm = TRUE) + 0.1
    y_position_mean <- max(df_subset[[mean_category]], na.rm = TRUE) + 0.1
    
    # Select color
    fill_color <- color_list[(variable_index - 1) %% length(color_list) + 1]
    
    # Data for count text
    count_df_var <- data.frame(
      group = names(group_counts),
      y = y_text_var,
      label = as.character(group_counts)
    )
    
    count_df_mean <- data.frame(
      group = names(group_counts),
      y = y_text_mean,
      label = as.character(group_counts)
    )
    
    # --- Variability Plot ---
    var_plot <- ggplot(df_subset, aes(x = !!sym(variable), y = log10(!!sym(var_category)))) +
      geom_boxplot(fill = fill_color, color = "darkgrey") +
      geom_text(data = count_df_var,
                aes(x = group, y = y, label = label),
                inherit.aes = FALSE, size = 3) +
      geom_signif(
        comparisons = list(c("No", "Yes")),
        annotations = sig_label,
        y_position = y_position_var,
        tip_length = 0.01,
        textsize = 2.5
      ) +
      ylab(paste0("log10(",condition," variability)")) +
      xlab(X_axis_name) +
      scale_x_discrete(labels = count_labels) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      theme_classic(base_size = 10) +
      theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)
      )
    
    # --- Mean Plot ---
    mean_plot <- ggplot(df_subset, aes(x = !!sym(variable), y = !!sym(mean_category))) +
      geom_boxplot(fill = fill_color, color = "darkgrey") +
      geom_text(data = count_df_mean,
                aes(x = group, y = y, label = label),
                inherit.aes = FALSE, size = 3) +
      geom_signif(
        comparisons = list(c("No", "Yes")),
        annotations = mean_sig_label,
        y_position = y_position_mean,
        tip_length = 0.01,
        textsize = 2.5
      ) +
      ylab(paste0("Mean transcript level (", condition, ")")) +
      xlab(X_axis_name) +
      scale_x_discrete(labels = count_labels) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      theme_classic(base_size = 10) +
      theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)
      )
    
    # Store plots
    varplot_list[[var_category]] <- var_plot
    meanplot_list[[mean_category]] <- mean_plot
  }
  
  # Save plots
  var_plots <- wrap_plots(varplot_list, nrow = 2)
  ggsave(plot = var_plots, filename = paste0("VarCor_", variable, "_plot.svg"), width = 2, height = 6, bg = "transparent",dpi=300)
  
  mean_plots <- wrap_plots(meanplot_list, nrow = 2)
  ggsave(plot = mean_plots, filename = paste0("MeanCor_", variable, "_plot.svg"), width = 2, height = 6, bg = "transparent",dpi=300)
  
  print(variable)
}

# View the correlation table
View(Categorical_test_table)

#################################################################
##### (3) Essentiality - Fraction plots and rank correlations ###
#################################################################

Essentiality_correlates <- c("essential","conditional_essential")

Y_axis_names <- c("% genes are always essential",
                  "% genes are conditionally essential")

Essentiality_correlation_table <- data.frame(
  Condition = character(),
  Correlate = character(),
  Number_of_genes = numeric(),
  Var_Spearman_Coefficient = numeric(),
  Var_Spearman_P_Value = numeric(),
  Mean_Spearman_Coefficient = numeric(),
  Mean_Spearman_P_value = numeric(),
  stringsAsFactors = FALSE
)

for (correlate_index in 1:length(Essentiality_correlates)){
  
  correlate <- Essentiality_correlates[correlate_index]
  meanplot_list <- list()
  varplot_list <- list()
  plot_list <- list()
  for (condition in c("Ctrl","HS")){
    
    var_category <- paste0(condition,"_",var_metric)
    mean_category <- paste0(condition,"_",mean_metric)
    
    df_subset <- Variability_correlates_df %>%
      dplyr::select(!!sym(var_category),!!sym(mean_category), !!sym(correlate)) %>%
      na.omit()
    
    df_subset[[var_category]] <- as.numeric(df_subset[[var_category]])
    df_subset[[mean_category]] <- as.numeric(df_subset[[mean_category]])
    
    # Create quantiles and percent essential per quantile
    var_quantile <- paste0(var_category," quantile")
    df_subset <- make_quantiles(data=df_subset,
                                quantile_column = df_subset[[var_category]],
                                n_quantiles = 10,
                                name_of_quantiles = var_quantile)
    percent_essential <- paste0("Percent ",correlate, "- variability")
    df_subset <- percent_per_quantile(data = df_subset,
                                                      column_to_percent_name = correlate,
                                                      column_with_quantiles = df_subset[[var_quantile]],
                                                      new_column_name = percent_essential)
    mean_quantile <- paste0(mean_category," quantile")
    df_subset <- make_quantiles(data=df_subset,
                                quantile_column = df_subset[[mean_category]],
                                n_quantiles = 10,
                                name_of_quantiles = mean_quantile)
    mean_percent_essential <- paste0("Percent ",correlate, "- mean")
    df_subset <- percent_per_quantile(data = df_subset,
                                      column_to_percent_name = correlate,
                                      column_with_quantiles = df_subset[[mean_quantile]],
                                      new_column_name = mean_percent_essential)
    
    # Extract Spearman correlation coefficient and p-value
    spearman_test <- cor.test(y=df_subset[[percent_essential]], x=as.numeric(df_subset[[var_quantile]]), method = "spearman")
    spearman_coef <- spearman_test$estimate
    spearman_p_value <- spearman_test$p.value
    
    mean_spearman_test <- cor.test(y=df_subset[[mean_percent_essential]], x=as.numeric(df_subset[[mean_quantile]]), method = "spearman")
    mean_spearman_coef <- mean_spearman_test$estimate
    mean_spearman_p_value <- mean_spearman_test$p.value
    
    # Add the results to the correlation table
    Essentiality_correlation_table <- rbind(Essentiality_correlation_table, data.frame(
      Condition = condition,
      Correlate = correlate,
      Number_of_genes = nrow(df_subset),
      Var_Spearman_Coefficient = spearman_coef,
      Var_Spearman_P_Value = spearman_p_value,
      Mean_Spearman_Coefficient = mean_spearman_coef,
      Mean_Spearman_P_value = mean_spearman_p_value,
      stringsAsFactors = FALSE
    ))
    
    # Define text
    Y_axis_name <- Y_axis_names[[correlate_index]]
    X_axis_name_mean <- paste0(condition,' mean quantile')
    X_axis_name_var <- paste0(condition,' variability quantile')
    number_of_genes <- nrow(df_subset)
    stat_test_result <- paste0("Spearman rho = ",format_statistic(spearman_coef),"\n",
                               " p = ",format_statistic(spearman_p_value),"\n","Number of genes = ",number_of_genes)
    mean_stat_test_result <- paste0("Spearman rho = ",format_statistic(mean_spearman_coef),"\n",
                                    "p = ",format_statistic(mean_spearman_p_value),"\n",'Number of genes = ',number_of_genes)
    
    # Plot
    var_plot <- ggplot(df_subset, aes(y = !!sym(percent_essential), x = !!sym(var_quantile))) +
      geom_boxplot() +
      geom_smooth(method = "lm", col = "blue")+
      ylab(Y_axis_name)+
      xlab(X_axis_name_var)+
      annotate("text", y = min(df_subset[,percent_essential]), x = 1, label = stat_test_result, size = 2.35, hjust = 0,vjust = 0)+
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            plot.title = element_text(size = 5))+
      theme_classic()  
    
    # Plot
    mean_plot <- ggplot(df_subset, aes(y = !!sym(percent_essential), x = !!sym(mean_quantile))) +
      geom_boxplot() +
      geom_smooth(method = "lm", col = "blue")+
      ylab(Y_axis_name)+
      xlab(X_axis_name_mean)+
      annotate("text", y = min(df_subset[,percent_essential]), x = 1, label = mean_stat_test_result, size = 2.35, hjust = 0,vjust = 0)+
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            plot.title = element_text(size = 5))+
      theme_classic()  

    plot_list[[mean_category]] <- mean_plot
    plot_list[[var_category]] <- var_plot
    meanplot_list[[mean_category]] <- mean_plot
    varplot_list[[var_category]] <- var_plot
  }
  # Save the plots
  mean_plots <- wrap_plots(meanplot_list, nrow = 2)  # Arrange in 2 columns
  print(mean_plots)
  ggsave(plot = mean_plots,filename = paste0("MeanCor_",correlate,"plot.svg"),width=3,height=6,dpi=300)
  var_plots <- wrap_plots(varplot_list, nrow = 2)  # Arrange in 2 columns
  print(var_plots)
  ggsave(plot = var_plots,filename = paste0("VarCor_",correlate,"plot.svg"),width=3,height=6,dpi=300)
  print(paste0(correlate_index,correlate))
}

View(Essentiality_correlation_table)

# Add these to the numeric correlation table
Numeric_correlation_table_complete <- bind_rows(Numeric_correlation_table, 
                                                Essentiality_correlation_table)
View(Numeric_correlation_table_complete)

######################################################################
##### (4) Plot the features correlation coefficients in one plot #####
######################################################################

# Add a new column for the correlate headings to include the number of genes
Numeric_correlation_table_complete <- Numeric_correlation_table_complete %>%
  mutate(
    Feature_description = paste0(gsub("_", " ", Correlate), " (", Number_of_genes," genes)")
  )
# Let us rename some features
Numeric_correlation_table_complete$Feature_description <- gsub("major shape ind", "promoter shape broadness", Numeric_correlation_table_complete$Feature_description)
Numeric_correlation_table_complete$Feature_description <- gsub("averagepi", "average nucleotide diversity (pi)", Numeric_correlation_table_complete$Feature_description)
Numeric_correlation_table_complete$Feature_description <- gsub("gc", "GC %", Numeric_correlation_table_complete$Feature_description)
Numeric_correlation_table_complete$Feature_description <- gsub("dNdS", "dN/dS", Numeric_correlation_table_complete$Feature_description)
Numeric_correlation_table_complete$Feature_description <- gsub("num dhs conditions.prox", "# tested conditions with DHS", Numeric_correlation_table_complete$Feature_description)
Numeric_correlation_table_complete$Feature_description <- gsub("number of", "#", Numeric_correlation_table_complete$Feature_description)
Numeric_correlation_table_complete$Feature_description <- gsub("ial", "iality", Numeric_correlation_table_complete$Feature_description)
Numeric_correlation_table_complete$Feature_description <- gsub("h2 h", "expression heritability", Numeric_correlation_table_complete$Feature_description)
View(Numeric_correlation_table_complete)

# Split by condition
Num_corr_table_Ctrl <- subset(Numeric_correlation_table_complete,Condition=='Ctrl')
Num_corr_table_HS <- subset(Numeric_correlation_table_complete,Condition=='HS')

# Create a column ranking feature importance 
Num_corr_table_Ctrl$importance <- rank(abs(Num_corr_table_Ctrl$Var_Spearman_Coefficient))
Num_corr_table_Ctrl$importance_var <- rank(abs(Num_corr_table_Ctrl$Var_Spearman_Coefficient))
Num_corr_table_Ctrl$importance_mean <- rank(abs(Num_corr_table_Ctrl$Mean_Spearman_Coefficient))
Num_corr_table_HS$importance <- rank(abs(Num_corr_table_HS$Var_Spearman_Coefficient))
Num_corr_table_HS$importance_var <- rank(abs(Num_corr_table_HS$Var_Spearman_Coefficient))
Num_corr_table_HS$importance_mean <- rank(abs(Num_corr_table_HS$Mean_Spearman_Coefficient))

# Long format to differentiate top feature for mean and variability
Num_corr_table_Ctrl <- Num_corr_table_Ctrl %>%
  pivot_longer(cols = c('importance_var','importance_mean'))
Num_corr_table_HS <- Num_corr_table_HS %>%
  pivot_longer(cols = c('importance_var','importance_mean'))

# Put the mean and variability correlations as a single column
Num_corr_table_Ctrl$Correlation <- ifelse(Num_corr_table_Ctrl$name=='importance_mean',Num_corr_table_Ctrl$Mean_Spearman_Coefficient,Num_corr_table_Ctrl$Var_Spearman_Coefficient)
Num_corr_table_HS$Correlation <- ifelse(Num_corr_table_HS$name=='importance_mean',Num_corr_table_HS$Mean_Spearman_Coefficient,Num_corr_table_HS$Var_Spearman_Coefficient)

# Ctrl plot
Num_corr_table_Ctrl$name <- factor(
  Num_corr_table_Ctrl$name,
  levels = c("importance_var", "importance_mean")
)
Num_corr_table_Ctrl <- Num_corr_table_Ctrl[order(Num_corr_table_Ctrl$name, decreasing = TRUE), ]
Ctrl_plot_table <- ggplot(
  Num_corr_table_Ctrl,
  aes(
    y = reorder(Feature_description, importance),
    x = value,
    color = name,
    fill = name,
    size = abs(Correlation),
    shape = factor(sign(Correlation))
  )
) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  xlab("Correlation rank (Ctrl)") +
  ylab("") +
  
  # Correct legend labels and values with named colors
  scale_color_manual(
    values = c("importance_var" = "#611BB8", "importance_mean" = "#4F8E4D"),
    name = "Type",
    labels = c("Variability", "Mean")
  ) +
  scale_fill_manual(
    values = c("importance_var" = "#611BB8", "importance_mean" = "#4F8E4D"),
    name = "Type",
    labels = c("Variability", "Mean"),
    guide = FALSE
  ) +
  
  scale_size_continuous(name = "Abs. cor", range = c(1, 8)) +
  scale_shape_manual(
    values = c("-1" = 25, "1" = 24),
    name = "Corr. sign",
    labels = c("-", "+")
  ) +
  
  guides(
    color = guide_legend(override.aes = list(size = 8), order = 1),
    shape = guide_legend(override.aes = list(size = 3, fill = "black"), order = 2)
  ) +
  
  coord_cartesian(clip = "off") +
  theme(
    # LEGEND ON RIGHT (default)
    legend.position = "right",
    legend.justification = "top",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9, color = 'black', hjust = 0.95, vjust = 0.2),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 1),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)  )
Ctrl_plot_table

# HS plot
Num_corr_table_HS$name <- factor(
  Num_corr_table_HS$name,
  levels = c("importance_var", "importance_mean")
)
Num_corr_table_HS <- Num_corr_table_HS[order(Num_corr_table_HS$name, decreasing = TRUE), ]
HS_plot_table <- ggplot(
  Num_corr_table_HS,
  aes(
    y = reorder(Feature_description, importance),
    x = value,
    color = name,
    fill = name,
    size = abs(Correlation),
    shape = factor(sign(Correlation))
  )
) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  xlab("Correlation rank (HS)") +
  ylab("") +
  
  # Correct legend labels and values with named colors
  scale_color_manual(
    values = c("importance_var" = "#611BB8", "importance_mean" = "#4F8E4D"),
    name = "Type",
    labels = c("Variability", "Mean")
  ) +
  scale_fill_manual(
    values = c("importance_var" = "#611BB8", "importance_mean" = "#4F8E4D"),
    name = "Type",
    labels = c("Variability", "Mean"),
    guide = FALSE
  ) +
  
  scale_size_continuous(name = "Abs. cor", range = c(1, 8)) +
  scale_shape_manual(
    values = c("-1" = 25, "1" = 24),
    name = "Corr. sign",
    labels = c("-", "+")
  ) +
  
  guides(
    color = guide_legend(override.aes = list(size = 8), order = 1),
    shape = guide_legend(override.aes = list(size = 3, fill = "black"), order = 2)
  ) +
  
  coord_cartesian(clip = "off") +
  theme(
    # LEGEND ON RIGHT (default)
    legend.position = "right",
    legend.justification = "top",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9, color = 'black', hjust = 0.95, vjust = 0.2),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 1),
    plot.margin = margin(t = 10, r = 0, b = 10, l = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)  )
HS_plot_table

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD\\")
ggsave(filename = 'MeanVarFeatureRankCtrlPlot.svg',Ctrl_plot_table,width = 6,height=4,dpi=300)
ggsave(filename = 'MeanVarFeatureRankHSPlot.svg',HS_plot_table,width = 6,height=4,dpi=300)

##############################
##### Save the workspace #####
##############################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD\\")

# Top variability correlates
columns_to_keep <- c('gene_ID','Ctrl_VST_MAD','HS_VST_MAD','Ctrl_VST_Mean','HS_VST_Mean',
                     Numeric_correlates,Essentiality_correlates,binary_vars)
Variability_correlates_df_top <- Variability_correlates_df[,columns_to_keep]
write.csv(Variability_correlates_df_top,"Variability_correlates_df.csv")

# Correlation tables
write.csv(Categorical_test_table,"Categorical_test_table.csv")
write.csv(Numeric_correlation_table_complete,"Numeric_correlation_table.csv")
write.csv(Essentiality_correlation_table,"Essentiality_correlation_table.csv")

save.image(file='Variability_Correlates.RData') 
