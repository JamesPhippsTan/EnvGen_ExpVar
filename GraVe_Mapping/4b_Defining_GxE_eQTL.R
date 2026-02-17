# GxE eQTL 
# Using the GEMMA gene-SNP pairs result
# to identify an appropriate GxE cutoff
# using tensorQTL results 

# Last Updated: 22/10/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls(all = T))

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
#load(file='4b_Defining_GxE_eQTL.RData')

#######################
##### Datasets ########
#######################

# Gene body locations obtained from droseqtl.org - these match the dataset though are not necessarily the latest annotations on flybase
Gene_body_locations <-fread('droseQTLorg_genes_positions.csv',header = T,select = c(1,2,3,6),
                            col.names = c("Gene_id","GENE_START","GENE_END","GENE_CHR"))
Gene_body_locations$phenotype_id_CHR <- gsub("X", "23", Gene_body_locations$phenotype_id_CHR)

# Load all tensorQTL all eQTL per condition
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
Ctrl_cis_eqtl_all <- read.table('Ctrl_cis_eqtl_result/Ctrl.cis_qtl.txt.gz',header = T, fill=TRUE)
HS_cis_eqtl_all <- read.table('HS_cis_eqtl_result/HS.cis_qtl.txt.gz',header = T, fill=TRUE)
Ctrl_trans_eqtl_all <- read.table('Ctrl_trans_eqtl_result/Ctrl.trans_qtl_pairs.txt.gz',header = T, fill=TRUE)
HS_trans_eqtl_all <- read.table('HS_trans_eqtl_result/HS.trans_qtl_pairs.txt.gz',header = T, fill=TRUE)

# Load sig eQTL at FDR < 0.05
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
Ctrl_cis_eqtl_sig <- read.table('Ctrl_cis_eqtl_sig.txt',header = T, fill=TRUE)
HS_cis_eqtl_sig <- read.table('HS_cis_eqtl_sig.txt',header = T, fill=TRUE)
Ctrl_trans_eqtl_sig <- read.table('Ctrl_trans_eqtl_sig.txt',header = T, fill=TRUE)
HS_trans_eqtl_sig <- read.table('HS_trans_eqtl_sig.txt',header = T, fill=TRUE)

###########################################################
##### (1) Defining tensorQTL GxE cutoff determination #####
###########################################################

# How to define a single GxE SNP if its significance was tested twice, once for each environment.
# One way is a gene-SNP pair significant in one condition based on a stringent cutoff
# and insignificant in another condition based on a looser cutoff 

# If the loose cutoff dataframe is too large to operate on
# Another way to think of the problem is defining 'non-GxE' or 'shared' gene-SNP
# that is significant in one condition based on a stringent cutoff 
# and significant in another condition based on a looser cutoff
# A GxE gene-SNP will hence be the number of FDR significant SNPs
# WITHOUT the shared gene-SNP

# Ultimately, this method will define two sets of GxE gene-SNP (one per condition), and one set of shared gene-SNP

# To get these - we can run a list intersect of the following
# (1) a dataframe of FDR < 0.05 gene-SNP pairs - these are the significant SNPs defined in the R codes
# (2) a dataframe of gene-SNP pairs with p-values < 0.05 (or some value allows the number of pairs to be closer to the true number of GxE pairs)

# First, we need to define FDR < 0.05 gene-SNP pairs 
Ctrl_cis_eqtl_sig$phenotype_id_SNP_pair <- paste0(Ctrl_cis_eqtl_sig$phenotype_id,"_",Ctrl_cis_eqtl_sig$variant_id)
HS_cis_eqtl_sig$phenotype_id_SNP_pair <- paste0(HS_cis_eqtl_sig$phenotype_id,"_",HS_cis_eqtl_sig$variant_id)
Ctrl_trans_eqtl_sig$phenotype_id_SNP_pair <- paste0(Ctrl_trans_eqtl_sig$phenotype_id,"_",Ctrl_trans_eqtl_sig$variant_id)
HS_trans_eqtl_sig$phenotype_id_SNP_pair <- paste0(HS_trans_eqtl_sig$phenotype_id,"_",HS_trans_eqtl_sig$variant_id)

# While we can grab the nominal pvalues from the original files, the original BH-corrected p-values are absent
# Let us correct them again here - this has been previously done in 'eQTL_prep_for_veQTL.R'

# Cis tests - BH method
n_cis_tests <- nrow(Ctrl_cis_eqtl_all)
Ctrl_cis_eqtl_all$efdr <- p.adjust(Ctrl_cis_eqtl_all$pval_nominal,method = 'BH',n = n_cis_tests)
HS_cis_eqtl_all$efdr <- p.adjust(HS_cis_eqtl_all$pval_nominal,method = 'BH',n = n_cis_tests)

# Trans tests - BH method
n_trans_tests <- 8763*383710 - n_cis_tests
Ctrl_trans_eqtl_all$efdr <- p.adjust(Ctrl_trans_eqtl_all$pval,method = 'BH',n = n_trans_tests)
HS_trans_eqtl_all$efdr <- p.adjust(HS_trans_eqtl_all$pval,method = 'BH',n = n_trans_tests)

######################################
##### (2a) cis-eQTL GxE ##############
######################################

# Initialize the output dataframe
row_labels <- c("cis_FDR<0.2", "cis_FDR<0.1", "cis_FDR<0.05",
                "trans_FDR<0.2", "trans_FDR<0.1", "trans_FDR<0.05"
                , "cisandtrans_FDR<0.2", "cisandtrans_FDR<0.1", "cisandtrans_FDR<0.05")

eQTL_GxE_table <- data.frame(
  SharedCutoff = row_labels,
  GxE_Ctrl_genes = integer(length(row_labels)),
  GxE_HS_genes = integer(length(row_labels)),
  Shared_genes = integer(length(row_labels)),
  GxE_Ctrl_SNPs = integer(length(row_labels)),
  GxE_HS_SNPs = integer(length(row_labels)),
  Shared_SNPs = integer(length(row_labels)),
  stringsAsFactors = FALSE
)


### --------- (1) CIS: FDR thresholds ---------
FDR_cutoffs <- c(0.2, 0.1, 0.05)
for (FDR in FDR_cutoffs) {
  row_name <- paste0("cis_FDR<", FDR)
  
  Ctrl_cis_cutoff <- subset(Ctrl_cis_eqtl_all, efdr < FDR)
  HS_cis_cutoff <- subset(HS_cis_eqtl_all, efdr < FDR)
  
  Shared_genes <- length(unique(c(
    intersect(Ctrl_cis_eqtl_sig$phenotype_id, HS_cis_cutoff$phenotype_id),
    intersect(HS_cis_eqtl_sig$phenotype_id, Ctrl_cis_cutoff$phenotype_id)
  )))
  
  GxE_Ctrl_genes <- length(setdiff(Ctrl_cis_eqtl_sig$phenotype_id, HS_cis_cutoff$phenotype_id))
  GxE_HS_genes   <- length(setdiff(HS_cis_eqtl_sig$phenotype_id, Ctrl_cis_cutoff$phenotype_id))
  
  Shared_SNPs <- length(unique(c(
    intersect(Ctrl_cis_eqtl_sig$variant_id, HS_cis_cutoff$variant_id),
    intersect(HS_cis_eqtl_sig$variant_id, Ctrl_cis_cutoff$variant_id)
  )))
  
  GxE_Ctrl_SNPs <- length(setdiff(Ctrl_cis_eqtl_sig$variant_id, HS_cis_cutoff$variant_id))
  GxE_HS_SNPs   <- length(setdiff(HS_cis_eqtl_sig$variant_id, Ctrl_cis_cutoff$variant_id))
  
  eQTL_GxE_table[eQTL_GxE_table$SharedCutoff == row_name, 2:7] <- c(
    GxE_Ctrl_genes, GxE_HS_genes, Shared_genes,
    GxE_Ctrl_SNPs, GxE_HS_SNPs, Shared_SNPs
  )
}


### --------- (2) TRANS: FDR thresholds ---------
for (FDR in FDR_cutoffs) {
  row_name <- paste0("trans_FDR<", FDR)

  Ctrl_trans_cutoff <- subset(Ctrl_trans_eqtl_all, efdr < FDR)
  HS_trans_cutoff   <- subset(HS_trans_eqtl_all, efdr < FDR)
  
  Ctrl_trans_cutoff_gene <- Ctrl_trans_cutoff$phenotype_id
  HS_trans_cutoff_gene   <- HS_trans_cutoff$phenotype_id
  Ctrl_trans_cutoff_snp  <- Ctrl_trans_cutoff$variant_id
  HS_trans_cutoff_snp    <- HS_trans_cutoff$variant_id
  
  Shared_genes <- length(unique(c(
    intersect(Ctrl_trans_eqtl_sig$phenotype_id, HS_trans_cutoff_gene),
    intersect(HS_trans_eqtl_sig$phenotype_id, Ctrl_trans_cutoff_gene)
  )))
  GxE_Ctrl_genes <- length(setdiff(Ctrl_trans_eqtl_sig$phenotype_id, HS_trans_cutoff_gene))
  GxE_HS_genes   <- length(setdiff(HS_trans_eqtl_sig$phenotype_id, Ctrl_trans_cutoff_gene))
  
  Shared_SNPs <- length(unique(c(
    intersect(Ctrl_trans_eqtl_sig$variant_id, HS_trans_cutoff_snp),
    intersect(HS_trans_eqtl_sig$variant_id, Ctrl_trans_cutoff_snp)
  )))
  GxE_Ctrl_SNPs <- length(setdiff(Ctrl_trans_eqtl_sig$variant_id, HS_trans_cutoff_snp))
  GxE_HS_SNPs   <- length(setdiff(HS_trans_eqtl_sig$variant_id, Ctrl_trans_cutoff_snp))
  
  eQTL_GxE_table[eQTL_GxE_table$SharedCutoff == row_name, 2:7] <- c(
    GxE_Ctrl_genes, GxE_HS_genes, Shared_genes,
    GxE_Ctrl_SNPs, GxE_HS_SNPs, Shared_SNPs
  )
}

### --------- (3) CISTRANS: FDR thresholds ---------
for (FDR in FDR_cutoffs) {
  row_name <- paste0("cisandtrans_FDR<", FDR)
  
  # New cutoff definitions
  Ctrl_cis_cutoff <- subset(Ctrl_cis_eqtl_all, efdr < FDR)
  HS_cis_cutoff <- subset(HS_cis_eqtl_all, efdr < FDR)
  Ctrl_trans_cutoff <- subset(Ctrl_trans_eqtl_all, efdr < FDR)
  HS_trans_cutoff   <- subset(HS_trans_eqtl_all, efdr < FDR)
  Ctrl_trans_cutoff_gene <- Ctrl_trans_cutoff$phenotype_id
  HS_trans_cutoff_gene   <- HS_trans_cutoff$phenotype_id
  Ctrl_trans_cutoff_snp  <- Ctrl_trans_cutoff$variant_id
  HS_trans_cutoff_snp    <- HS_trans_cutoff$variant_id
  
  # Get the intersects
  Shared_genes <- length(unique(c(
    intersect(c(Ctrl_cis_eqtl_sig$phenotype_id,Ctrl_trans_eqtl_sig$phenotype_id), 
              c(HS_cis_cutoff$phenotype_id,HS_trans_cutoff_gene)),
    intersect(c(HS_cis_eqtl_sig$phenotype_id,HS_trans_eqtl_sig$phenotype_id), 
              c(Ctrl_cis_cutoff$phenotype_id,Ctrl_trans_cutoff_gene))
  )))
  
  GxE_Ctrl_genes <- length(setdiff(c(Ctrl_cis_eqtl_sig$phenotype_id,Ctrl_trans_eqtl_sig$phenotype_id), 
                                   c(HS_cis_cutoff$phenotype_id,HS_trans_cutoff_gene)))
  GxE_HS_genes   <- length(setdiff(c(HS_cis_eqtl_sig$phenotype_id,HS_trans_eqtl_sig$phenotype_id), 
                                   c(Ctrl_cis_cutoff$phenotype_id,Ctrl_trans_cutoff_gene)))
  
  Shared_SNPs <- length(unique(c(
    intersect(c(Ctrl_cis_eqtl_sig$variant_id,Ctrl_trans_eqtl_sig$variant_id), 
              c(HS_cis_cutoff$variant_id,HS_trans_cutoff_snp)),
    intersect(c(HS_cis_eqtl_sig$variant_id,HS_trans_eqtl_sig$variant_id), 
              c(Ctrl_cis_cutoff$variant_id,Ctrl_trans_cutoff_snp))
  )))
  
  GxE_Ctrl_SNPs <- length(setdiff(c(Ctrl_cis_eqtl_sig$variant_id,Ctrl_trans_eqtl_sig$variant_id), 
                                  c(HS_cis_cutoff$variant_id,HS_trans_cutoff_snp)))
  GxE_HS_SNPs   <- length(setdiff(c(HS_cis_eqtl_sig$variant_id,HS_trans_eqtl_sig$variant_id), 
                                  c(Ctrl_cis_cutoff$variant_id,Ctrl_trans_cutoff_snp)))
  
  eQTL_GxE_table[eQTL_GxE_table$SharedCutoff == row_name, 2:7] <- c(
    GxE_Ctrl_genes, GxE_HS_genes, Shared_genes,
    GxE_Ctrl_SNPs, GxE_HS_SNPs, Shared_SNPs
  )
}


### View Final Summary Table for eQTL
View(eQTL_GxE_table)

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\GxE")
write.csv(eQTL_GxE_table,'eQTL_GxE_table.csv',row.names = F)

# Turn it into a percentage
eQTL_GxE_table_percent <- eQTL_GxE_table
eQTL_GxE_table_percent[1:3,2:4] <- eQTL_GxE_table[1:3,2:4]*100/(sum(eQTL_GxE_table[1,2:4]))
eQTL_GxE_table_percent[4:6,2:4] <- eQTL_GxE_table[4:6,2:4]*100/(sum(eQTL_GxE_table[5,2:4]))
eQTL_GxE_table_percent[7:9,2:4] <- eQTL_GxE_table[7:9,2:4]*100/(sum(eQTL_GxE_table[9,2:4]))
eQTL_GxE_table_percent[1:3,5:7] <- eQTL_GxE_table[1:3,5:7]*100/(sum(eQTL_GxE_table[1,5:7]))
eQTL_GxE_table_percent[4:6,5:7] <- eQTL_GxE_table[4:6,5:7]*100/(sum(eQTL_GxE_table[5,5:7]))
eQTL_GxE_table_percent[7:9,5:7] <- eQTL_GxE_table[7:9,5:7]*100/(sum(eQTL_GxE_table[9,5:7]))
eQTL_GxE_table_percent[,2:7] <- round(eQTL_GxE_table_percent[,2:7],digits = 2)
View(eQTL_GxE_table_percent)
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\GxE")
write.csv(eQTL_GxE_table_percent,'eQTL_GxE_table_percent.csv')

# Make a plot of the percentages
eQTL_GxE_table_percent_separated <- eQTL_GxE_table_percent %>%
  separate(SharedCutoff, into = c("cis_or_trans", "shared_cutoff"), sep = "_")
eQTL_GxE_table_percent_separated$shared_cutoff <- factor(
  eQTL_GxE_table_percent_separated$shared_cutoff,
  levels = c("FDR<0.05", "FDR<0.1", "FDR<0.2")
)

# Order factors
eQTL_GxE_table_percent_separated$shared_cutoff <- factor(eQTL_GxE_table_percent_separated$shared_cutoff,
                                                               levels = c("FDR<0.05", "FDR<0.1", "FDR<0.2")
                                                        )

# Define the columns to plot
columns_to_plot <- colnames(eQTL_GxE_table_percent_separated)[3:8]


# Loop over columns
for (col_name in columns_to_plot) {
  p <- ggplot(eQTL_GxE_table_percent_separated,
              aes(x = shared_cutoff, y = .data[[col_name]],
                  color = cis_or_trans, group = cis_or_trans)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    theme_classic() +
    labs(
      title = 'eQTL',
      x = "Significance cutoff in other condition",
      y = gsub('_',' ',paste0(col_name, " (%)")),
      color = "Regulation Type"
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot
  ggsave(filename = paste0(col_name, "_eQTL_dotplot.svg"),
         plot = p, width = 6, height = 4, dpi = 300)
}

##############################################################
##### Save files and working environment so far ##############
##############################################################

save.image(file='4b_Defining_GxE_eQTL.RData')

