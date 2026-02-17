# Investigating whether eQTL or veQTL are enriched in functional genomic regions

# Last Updated: 10/9/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls(all = T))

library(tibble)
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Regulatory_Feature_Enrichments")
load(file = '7b_Regulatory_Feature_Enrichments.RData')

# Functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('gProfiler_Functions.R')
source('Quantile_Functions.R')
source('Variability_Functions.R')
source('QTL_Analysis_Functions.R')

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code\\GraVe_Mapping")

#######################
##### Datasets ########
#######################

# Load all SNPs
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
all_SNPs <- read.table('SNPs_original_positions.txt')

# Load cis-QTL SNPs (for downstream testing)
cis_QTL_df <- read.table("HS_cis_veQTL.txt", header = T,skipNul = T)
all_cis_SNPs <- unique(paste0(sub(23,"X",cis_QTL_df$CHROM),cis_QTL_df$POS))

# Load hotspot SNPs and summary for reference
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs")
hotspot_SNPs <- read.csv('hotspot_SNPs.csv')
hotspot_SNPs <- hotspot_SNPs[,-1]
hotspot_SNPs_summaries <- read.csv('hotspot_SNPs_summaries.csv')
hotspot_SNPs_summaries_short <- hotspot_SNPs_summaries[,c('X','Total_number')]
hotspot_SNPs_summaries_short$X <- gsub('vGenes','veQTL',hotspot_SNPs_summaries_short$X) 
hotspot_SNPs_summaries_short$X <- gsub('eGenes','eQTL',hotspot_SNPs_summaries_short$X) 
View(hotspot_SNPs_summaries_short)
hotspot_SNPs_summary <- data.frame(row.names = hotspot_SNPs_summaries_short$X,
                                   Total_number=hotspot_SNPs_summaries_short$Total_number)

# Significant eQTL 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
Ctrl_cis_eQTL <-read.table("Ctrl_cis_eqtl_sig.txt",header = T)
Ctrl_cis_eQTL$SNP <- Ctrl_cis_eQTL$variant_id
HS_cis_eQTL <-read.table("HS_cis_eqtl_sig.txt",header = T)
HS_cis_eQTL$SNP <- HS_cis_eQTL$variant_id
Ctrl_trans_eQTL <-read.table("Ctrl_trans_eqtl_sig.txt",header = T)
Ctrl_trans_eQTL$SNP <- Ctrl_trans_eQTL$variant_id
HS_trans_eQTL <-read.table("HS_trans_eqtl_sig.txt",header = T)
HS_trans_eQTL$SNP <- HS_trans_eQTL$variant_id

# Significant veQTL
Ctrl_cis_veQTL <-read.csv("Ctrl_cis_veQTL_sig.csv")
HS_cis_veQTL <-read.csv("HS_cis_veQTL_sig.csv")
Ctrl_trans_veQTL <-read.csv("Ctrl_trans_veQTL_sig.csv")
HS_trans_veQTL <-read.csv("HS_trans_veQTL_sig.csv")

# Split the SNP names into SNP chromosomes and positions
Ctrl_cis_eQTL <-extract_chromosome_and_pos(Ctrl_cis_eQTL,input_col = 'SNP',keep_original = T)
HS_cis_eQTL <-extract_chromosome_and_pos(HS_cis_eQTL,input_col = 'SNP',keep_original = T)
Ctrl_trans_eQTL <-extract_chromosome_and_pos(Ctrl_trans_eQTL,input_col = 'SNP',keep_original = T)
HS_trans_eQTL <-extract_chromosome_and_pos(HS_trans_eQTL,input_col = 'SNP',keep_original = T)
Ctrl_cis_veQTL <-extract_chromosome_and_pos(Ctrl_cis_veQTL,input_col = 'SNP',keep_original = T)
HS_cis_veQTL <-extract_chromosome_and_pos(HS_cis_veQTL,input_col = 'SNP',keep_original = T)
Ctrl_trans_veQTL <-extract_chromosome_and_pos(Ctrl_trans_veQTL,input_col = 'SNP',keep_original = T)
HS_trans_veQTL <-extract_chromosome_and_pos(HS_trans_veQTL,input_col = 'SNP',keep_original = T)

# Load SNP effects - introns, exons, etc. 
# These already have a per-SNP basis
SNP_effects <- read.csv("SnpEff_Annotations_final.csv",header = T,row.names = 1)
SNP_effects <- SNP_effects[SNPs_position_map$SNP,]
SNP_effects <- rownames_to_column(SNP_effects,var = "SNP")
names(SNP_effects) <- gsub("^X", "", names(SNP_effects))

# Promoters and coordinates
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\EukaryoticPromoterDb")
Promoters <-read.table("drosophila_epdnew_HcTKY.bed")
colnames(Promoters) <- c('CHROM','START','POS','gene_name','1','STRAND')
Promoters$CHROM = str_remove(Promoters$CHROM, "^chr")

# scEnhancers and coordinates
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\scEnhancers\\scEnhancer_and_GeneExpression")
scEnhancers <- process_scEnhancer_tables(getwd())
scEnhancers$element_name <- scEnhancers$name

# Transcription factor binding sites (TFBSs) and coordinates
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\REDfly_Data")
TFBSs <-read_tsv("REDfly_tfbs.tsv")
TFBSs$Chromosome = str_remove(TFBSs$Chromosome, "^chr")
colnames(TFBSs)
Converted_Colnames <- c('CHROM','START','END','element_name','gene_identifier')
Col_to_convert_names <- c(1,2,3,4,8)
colnames(TFBSs)[Col_to_convert_names] <- Converted_Colnames
colnames(TFBSs)

# Cis-regulatory modules (CRMs - e.g., enhancers and silencers) and coordinates
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\REDfly_Data")
CRMs <-read_tsv("REDfly_crm.tsv")
CRMs$Chromosome = str_remove(CRMs$Chromosome, "^chr")
Converted_Colnames <- c('CHROM','START','END','element_name','gene_identifier')
Col_to_convert_names <- c(1,2,3,4,7)
colnames(CRMs)
colnames(CRMs)[Col_to_convert_names] <- Converted_Colnames
colnames(CRMs)

# Create subsets of CRMs
CRMs_enhancers <- CRMs[grepl("enhancer", CRMs$`Enhancer/Silencer Status`), ]
CRMs_silencers <- CRMs[grepl("silencer", CRMs$`Enhancer/Silencer Status`), ]
CRMs_others <- CRMs[is.na(CRMs$`Enhancer/Silencer Status`), ]

CRMs_female <- CRMs[grepl("f|both", CRMs$Sex), ]
CRMs_female_enhancers <-  CRMs_female[grepl("enhancer", CRMs_female$`Enhancer/Silencer Status`), ]
CRMs_female_silencers <- CRMs_female[grepl("silencer", CRMs_female$`Enhancer/Silencer Status`), ]
# Note - no others and mostly enhancers

CRMs_female_adult <- CRMs_female[grepl("adult", CRMs_female$`Stage On`), ]
CRMs_female_adult_enhancers <-  CRMs_female_adult[grepl("enhancer", CRMs_female_adult$`Enhancer/Silencer Status`), ]
CRMs_female_adult_silencers <- CRMs_female_adult[grepl("silencer", CRMs_female_adult$`Enhancer/Silencer Status`), ]

#############################################
##### (1) Test if the functions work ########
#############################################

# Test on expected positives

# Search the list of promoters
head(Promoters)
known_TSSs <- c('443759','470499')
promoter_search(snp_list = known_TSSs,
                promoter_df = Promoters,
                window_size = 50)

# Search the list of scEnhancers
head(scEnhancers)
known_enhancers <- c('2L9963644','2L9966900')
genomic_element_search(snp_list = known_enhancers,genomic_element_df = scEnhancers)

# Search the list of TFBSs
head(TFBSs)
known_TFBSs <- c('2R11264951','3L5813857')
genomic_element_search(snp_list = known_TFBSs,genomic_element_df = TFBSs)

# Search the list of CRMs - including enhancers and silencers
head(CRMs)
known_CRMs <- c('2R14801005','2R16639869')
genomic_element_search(snp_list = known_CRMs, genomic_element_df = CRMs)
genomic_element_search(snp_list = known_CRMs, genomic_element_df = CRMs_enhancers)
genomic_element_search(snp_list = known_CRMs, genomic_element_df = CRMs_female)
head(CRMs_female_adult)
known_female_adult_CRMs <- c('3L17822836','2L11929570')
genomic_element_search(snp_list = known_female_adult_CRMs, genomic_element_df = CRMs_female_adult)
genomic_element_search(snp_list = known_female_adult_CRMs, genomic_element_df = CRMs_female_adult_enhancers)

#############################################
##### (2) Assign a feature to each SNP ######
#############################################

# All SNPs
all_SNP_IDs <- all_SNPs$V3

# Search the list of TFBSs
#all_SNP_TFBSs <- genomic_element_search(snp_list = all_SNP_IDs,genomic_element_df = TFBSs)

# Search the list of scEnhancers
#all_SNP_scEnhancers <- genomic_element_search(snp_list = all_SNP_IDs,genomic_element_df = scEnhancers)

# Search the list of promoters - 50 bp window
#all_SNP_promoters_50 <- promoter_search(snp_list = all_SNP_IDs,
#                                        promoter_df = Promoters,
#                                        window_size = 50)

# Search the list of promoters - 100 bp window
#all_SNP_promoters_100 <- promoter_search(snp_list = all_SNP_IDs,
#                                         promoter_df = Promoters,
#                                         window_size = 100)

# Search the list of CRMs
#all_SNP_CRMs <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs)

# CRM subsets
#all_SNP_CRMs_enhancers <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs_enhancers)
#all_SNP_CRMs_silencers <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs_silencers)
#all_SNP_CRMs_others <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs_others)
#all_SNP_CRMs_female <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs_female)
#all_SNP_CRMs_female_enhancers <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs_female_enhancers)
#all_SNP_CRMs_female_silencers <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs_female_silencers)
#all_SNP_CRMs_female_adult <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs_female_adult)
#all_SNP_CRMs_female_adult_enhancers <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs_female_adult_enhancers)
#all_SNP_CRMs_female_adult_silencers <- genomic_element_search(snp_list = all_SNP_IDs, genomic_element_df = CRMs_female_adult_silencers)

# Collapsed versions of the TFBS and CRM dataframes
collapse_df <- function(searched_df){
  searched_df %>%
    group_by(snp_name) %>%
    summarise(
      element_name = paste(unique(element_name), collapse = "; "),
      gene_identifier = paste(unique(gene_identifier), collapse = "; "),
      overlaps = paste(unique(overlaps), collapse = "; "),
      .groups = "drop")
  }

all_SNP_TFBSs_collapsed <- collapse_df(all_SNP_TFBSs)
all_SNP_CRMs_collapsed <- collapse_df(all_SNP_CRMs)
all_SNP_CRMs_enhancers_collapsed <- collapse_df(all_SNP_CRMs_enhancers)
all_SNP_CRMs_silencers_collapsed <- collapse_df(all_SNP_CRMs_silencers) 
all_SNP_CRMs_others_collapsed <- collapse_df(all_SNP_CRMs_others) 
all_SNP_CRMs_female_collapsed <- collapse_df(all_SNP_CRMs_female) 
all_SNP_CRMs_female_enhancers_collapsed <- collapse_df(all_SNP_CRMs_female_enhancers) 
all_SNP_CRMs_female_silencers_collapsed <- collapse_df(all_SNP_CRMs_female_silencers) 
all_SNP_CRMs_female_adult_collapsed <- collapse_df(all_SNP_CRMs_female_adult) 
all_SNP_CRMs_female_adult_enhancers_collapsed <- collapse_df(all_SNP_CRMs_female_adult_enhancers) 
all_SNP_CRMs_female_adult_silencers_collapsed <- collapse_df(all_SNP_CRMs_female_adult_silencers) 

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Hotspot_SNPs")


###############################################################################
#### (3) Do top eSNPs and vSNPs tend to fall in genomic functional regions ####
###############################################################################

# Create single dataframe of all assayed SNPs and their features
all_SNP_features_reg <- data.frame(SNP = all_SNP_IDs,
                          Promoter = all_SNP_promoters_100$near_tss,
                          CRMs = all_SNP_CRMs_collapsed$overlaps,
                          Silencers = all_SNP_CRMs_silencers_collapsed$overlaps,
                          Enhancers = all_SNP_CRMs_enhancers_collapsed$overlaps,
                          CRMs_Other = all_SNP_CRMs_others_collapsed$overlaps,
                          CRMs_female_adult = all_SNP_CRMs_female_adult_collapsed$overlaps,
                          Enhancers_female_adult = all_SNP_CRMs_female_adult_enhancers_collapsed$overlaps,
                          Silencers_female_adult = all_SNP_CRMs_female_adult_silencers_collapsed$overlaps,
                          TFBS = all_SNP_TFBSs_collapsed$overlaps,
                          scEnhancers = all_SNP_scEnhancers$overlaps)
all_SNP_features <- merge(all_SNP_features_reg,SNP_effects,by='SNP')

# Initialize a list to store results
results <- list()

# Get the features to test for enrichment
SNP_features <- colnames(all_SNP_features)[2:ncol(all_SNP_features)]

# All cis- and trans-SNPs
# Used to create the background SNPs for the chi-square testing
all_trans_SNPs_features <- column_to_rownames(all_SNP_features,var = 'SNP')
all_cis_SNPs_features <- all_trans_SNPs_features[all_cis_SNPs,] 

# Initialize a list to store results
results <- list()

# List of condition-specific cis-SNP data frames
cis_condition_SNPs <- list(
  Ctrl_cis_eQTL = unique(Ctrl_cis_eQTL$SNP),
  HS_cis_eQTL = unique(HS_cis_eQTL$SNP),
  Ctrl_cis_veQTL= unique(Ctrl_cis_veQTL$SNP),
  HS_cis_veQTL = unique(HS_cis_veQTL$SNP)
)

# List of condition-specific trans-SNP data frames
trans_condition_SNPs <- list(
  Ctrl_trans_eQTL = unique(Ctrl_trans_eQTL$SNP),
  HS_trans_eQTL = unique(HS_trans_eQTL$SNP),
  Ctrl_trans_veQTL = unique(Ctrl_trans_veQTL$SNP),
  HS_trans_veQTL = unique(HS_trans_veQTL$SNP)
)

# Loop through each column
for (type_col in SNP_features) {
  
  # Perform tests for cis-SNPs
  for (cond_name in names(cis_condition_SNPs)) {
    condition_SNPs <- cis_condition_SNPs[[cond_name]]
    result  <- create_contingency_table(condition_SNPs, all_cis_SNPs_features, type_col)
    test_result <- perform_chi_square_test(result$contingency_table)
    
    # Store the result
    results[[paste0(cond_name, "_", type_col)]] <- list(
      cond = cond_name,
      type = type_col,
      p_value = test_result$p.value,
      condition_ratio = result$condition_ratio,
      remaining_ratio = result$remaining_ratio,
      condition_counts = result$condition_counts,
      remaining_counts = result$remaining_counts   
    )
  }
  
  # Perform tests for trans-SNPs
  for (cond_name in names(trans_condition_SNPs)) {
    condition_SNPs <- trans_condition_SNPs[[cond_name]]
    result <- create_contingency_table(condition_SNPs, all_trans_SNPs_features, type_col)
    test_result <- perform_chi_square_test(result$contingency_table)
    
    # Store the result
    results[[paste0(cond_name, "_", type_col)]] <- list(
      cond = cond_name,
      type = type_col,
      p_value = test_result$p.value,
      condition_ratio = result$condition_ratio,
      remaining_ratio = result$remaining_ratio,
      condition_counts = result$condition_counts,
      remaining_counts = result$remaining_counts   
    )  }
}

# Create a summary table
summary_table <- data.frame(
  Condition = sapply(results, function(x) x$cond),
  Feature = sapply(results, function(x) x$type),
  P_Value = sapply(results, function(x) x$p_value),
  Subset_Ratio = sapply(results, function(x) x$condition_ratio),
  Background_Ratio = sapply(results, function(x) x$remaining_ratio),
  ObsExp = sapply(results, function(x) x$condition_ratio/x$remaining_ratio),
  Cond_yes_fraction = sapply(results, function(x) x$condition_counts[["yes"]]/(x$condition_counts[["no"]]+x$condition_counts[["yes"]]))
)

# Apply multiple testing correction and format in a way that is visible
summary_table$BH_pval <- p.adjust(summary_table$P_Value,method = 'BH')
summary_table$pval_short <- sapply(summary_table$P_Value, function(x) format_statistic(as.numeric(x)))
summary_table$BH_pval_short <- sapply(summary_table$BH_pval, function(x) format_statistic(as.numeric(x)))
View(summary_table)

# Add a column for the direction of the difference
summary_table$Direction <- ifelse(
  summary_table$BH_pval < 0.05,
  ifelse(
    summary_table$Subset_Ratio > summary_table$Background_Ratio,
    paste0("More (p=",summary_table$BH_pval_short,")"),
    paste0("Fewer (p=",summary_table$BH_pval_short,")")
  ),
  paste0("N.D. (p=",summary_table$BH_pval_short,")")
)

# Order to make the eQTL-veQTL comparisons
Category_order = c('Ctrl_cis_eQTL', 
                   'HS_cis_eQTL',
                   'Ctrl_trans_eQTL', 
                   'HS_trans_eQTL',
                   'Ctrl_cis_veQTL', 
                   'HS_cis_veQTL',
                   'Ctrl_trans_veQTL', 
                   'HS_trans_veQTL') 

# Make a table on the direction of change and significance of the chi square test
Feature_enrichment_table <- summary_table[,c("Condition","Feature","Direction")] %>%
  pivot_wider(names_from = Condition, values_from = Direction)
Feature_enrichment_table <- column_to_rownames(as.data.frame(Feature_enrichment_table),var = 'Feature')
Feature_enrichment_table_T <- t(Feature_enrichment_table[,Category_order])
Feature_enrichment_table_T <- cbind(Feature_enrichment_table_T,
                                                Total_SNPs = hotspot_SNPs_summary[Category_order,'Total_number'])
View(Feature_enrichment_table_T)

# More detailed version of the above, including the yes:no fraction of the qtl vs the non-qtls
Feature_enrichment_table_ObsExpRatio<- summary_table[,c("Condition","Feature","ObsExp")] %>%
  pivot_wider(names_from = Condition, values_from = ObsExp)
Feature_enrichment_table_ObsExpRatio<- column_to_rownames(as.data.frame(Feature_enrichment_table_ObsExpRatio),var = 'Feature')
Feature_enrichment_table_T_ObsExpRatio<- t(Feature_enrichment_table_ObsExpRatio[,Category_order])
Feature_enrichment_table_T_ObsExpRatio <- cbind(Feature_enrichment_table_T_ObsExpRatio,
                                                Total_SNPs = hotspot_SNPs_summary[Category_order,'Total_number'])
View(Feature_enrichment_table_T_ObsExpRatio)

# Including the yes:no fraction of the qtl 
Feature_enrichment_table_YesFraction <- summary_table[,c("Condition","Feature","Cond_yes_fraction")] %>%
  pivot_wider(names_from = Condition, values_from = Cond_yes_fraction)
Feature_enrichment_table_YesFraction <- column_to_rownames(as.data.frame(Feature_enrichment_table_YesFraction),var = 'Feature')
Feature_enrichment_table_T_YesFraction <- t(Feature_enrichment_table_YesFraction[,Category_order])
Feature_enrichment_table_T_YesFraction <- cbind(Feature_enrichment_table_T_YesFraction,
                                                Total_SNPs = hotspot_SNPs_summary[Category_order,'Total_number'])
View(Feature_enrichment_table_T_YesFraction)

##############################################
##### Save the working environment ###########
##############################################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Regulatory_Feature_Enrichments")

write.csv(x = Feature_enrichment_table_T,file = 'Regulatory_feature_enrichment_table.csv')
write.csv(x = Feature_enrichment_table_T_ObsExpRatio,file = 'ObsExpRatio_regulatory_feature_enrichment_table.csv')
write.csv(x = Feature_enrichment_table_T_YesFraction,file = 'YesFraction_regulatory_feature_enrichment_table.csv')

save.image(file='7_Regulatory_Feature_Enrichments.RData')

