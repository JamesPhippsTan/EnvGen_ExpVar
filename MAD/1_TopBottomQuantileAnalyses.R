# Transcriptome variability analyses of Adult Netherlands D.mel heads on normal vs high sugar diets
# Looking across quantiles

# The top and bottom quantiles may offer a biased view of the transcriptome

# Last Updated: 4/8/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(ggplot2)
library(ggvenn)
library(gridExtra)
library(grid)
library(data.table)

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD")
load(file='TopBottomQuantileAnalyses_MAD.RData')

# Functions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Code")
source('Quantile_Functions.R')
source('gProfiler_Functions.R')
source('Variability_Functions.R')

######################
##### Datasets #######
######################

# Dataset
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\")
MeanVarTestResults <- read.csv('Gene_MeanVar_Table.csv',row.names = 1) 

# Set background genes
background_genes <- rownames(MeanVarTestResults)
number_of_expressed_genes = length(background_genes)

####################################
##### 1) Create quantiles ##########
####################################

# Quantiles
MeanVarTestResults <- make_quantiles(MeanVarTestResults,
                                     quantile_column = MeanVarTestResults$Ctrl_VST_MAD,
                                     n_quantiles = 20,
                                     name_of_quantiles = 'VarCtrlQuantile')
MeanVarTestResults <- make_quantiles(MeanVarTestResults,
                                     quantile_column = MeanVarTestResults$HS_VST_MAD,
                                     n_quantiles = 20,
                                     name_of_quantiles = 'VarHSQuantile')
MeanVarTestResults <- make_quantiles(MeanVarTestResults,
                                     quantile_column = MeanVarTestResults$Ctrl_VST_Mean,
                                     n_quantiles = 20,
                                     name_of_quantiles = 'MeanCtrlQuantile')
MeanVarTestResults <- make_quantiles(MeanVarTestResults,
                                     quantile_column = MeanVarTestResults$HS_VST_Mean,
                                     n_quantiles = 20,
                                     name_of_quantiles = 'MeanHSQuantile')

##################################################################
##### 2) Overlap of gene membership in the same quantile #########
##################################################################

# Look at the whole transcriptome correlation on gene-wise varaibility
Variability_by_Diet_Correlation <- cor.test(MeanVarTestResults$HS_VST_MAD,MeanVarTestResults$Ctrl_VST_MAD)
cor = format_statistic(Variability_by_Diet_Correlation$estimate)
pval = Variability_by_Diet_Correlation$p.value

Variability_by_Condition_Plot <- ggplot(MeanVarTestResults,aes(y=HS_VST_MAD,x=Ctrl_VST_MAD))+
  geom_point(color='darkgrey',alpha=0.8)+
  ylab('HS transcript level variability (MAD)')+
  xlab('Ctrl transcript level variability (MAD)')+
  annotate("text", x = 0, y = 1.05, label = paste0("Pearson correlation = ",cor,"\n","p-value < 2.2e-16"), size = 3.5, hjust = 0)+
  geom_smooth(method = "lm", col = "blue",linetype='dashed')+
  theme_classic()
Variability_by_Condition_Plot
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD")
ggsave(dpi = 300,plot = Variability_by_Condition_Plot,filename = "Variability_by_Condition_Plot.svg",width = 3,height = 3)

Mean_by_Diet_Correlation <- cor.test(MeanVarTestResults$HS_VST_Mean,MeanVarTestResults$Ctrl_VST_Mean)
cor_mean = format_statistic(Mean_by_Diet_Correlation$estimate)
pval_mean = Mean_by_Diet_Correlation$p.value

Mean_by_Condition_Plot <- ggplot(MeanVarTestResults,aes(y=HS_VST_Mean,x=Ctrl_VST_Mean))+
  geom_point(color='darkgrey',alpha=0.8)+
  ylab('HS mean transcript level')+
  xlab('Ctrl mean transcript level')+
  annotate("text", x = 4.5, y = Inf, label = paste0("Pearson correlation = ",cor_mean,"\n","p-value < 2.2e-16"), size = 3.5, hjust = 0,vjust=1)+
  geom_smooth(method = "lm", col = "blue",linetype='dashed')+
  theme_classic()
Mean_by_Condition_Plot
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Publication\\ggplots/")
ggsave(dpi = 300,plot = Mean_by_Condition_Plot,filename = "Mean_by_Condition_Plot.svg",width = 3,height = 3)


# Overlap = intersection/size of quantile (as opposed to size of whole gene set)
Quantile_overlap <- data.frame(Quantile = 1:20, 
                               Number_of_genes = NA*20, 
                               Var_Quant_Overlap = NA*20,
                               Mean_Quant_Overlap = NA*20)
for (quantile in 1:20){
  CtrlVarQ <- subset(MeanVarTestResults, VarCtrlQuantile == quantile)
  HSVarQ <- subset(MeanVarTestResults, VarHSQuantile == quantile)
  CtrlMeanQ <- subset(MeanVarTestResults, MeanCtrlQuantile == quantile)
  HSMeanQ <- subset(MeanVarTestResults, MeanHSQuantile == quantile)
  Quantile_overlap[quantile,'Number_of_genes'] <- nrow(HSVarQ)
  Quantile_overlap[quantile,'Var_Quant_Overlap'] <- 100*length(intersect(rownames(CtrlVarQ),rownames(HSVarQ)))/nrow(CtrlVarQ)
  Quantile_overlap[quantile,'Mean_Quant_Overlap'] <- 100*length(intersect(rownames(CtrlMeanQ),rownames(HSMeanQ)))/nrow(CtrlMeanQ)
}

Mean_Quant_Overlap_Plot <- ggplot(Quantile_overlap,aes(y=Mean_Quant_Overlap,x=Quantile))+
  geom_point(color='darkgrey')+
  ylab('% genes overlap between diets')+
  xlab('Mean quantile') + ylim(0,100)+
  theme_classic()
Mean_Quant_Overlap_Plot
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD")
ggsave(dpi = 300,plot = Mean_Quant_Overlap_Plot,
       filename = "Mean_Quant_Overlap_Plot.svg",
       height = 3,
       width=3)

Var_Quant_Overlap_Plot <- ggplot(Quantile_overlap,aes(y=Var_Quant_Overlap,x=Quantile))+
  geom_point(color='darkgrey')+
  ylab('% genes overlap between diets')+
  xlab('Variability quantile') + ylim(0,100)+
  theme_classic()
Var_Quant_Overlap_Plot
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD")
ggsave(dpi = 300,plot =Var_Quant_Overlap_Plot,
       filename ="Var_Quant_Overlap_Plot.svg",
       height = 3,
       width=3)

###########################################################################
##### 3) GO terms across the top and bottom variability quantiles #########
###########################################################################

# We are interested in genes that belong to the top and bottom quantile in both conditions
shared_or_different_quantile <- list()
quantiles <- c(1,20)
for (quantile in quantiles){
  CtrlVarQ <- subset(MeanVarTestResults, VarCtrlQuantile == quantile)
  HSVarQ <- subset(MeanVarTestResults, VarHSQuantile == quantile)
  shared_or_different_quantile[[paste0("VarCtrlQuantile",quantile)]] <- setdiff(rownames(CtrlVarQ),rownames(HSVarQ))
  shared_or_different_quantile[[paste0("VarHSQuantile",quantile)]] <- setdiff(rownames(HSVarQ),rownames(CtrlVarQ))
  shared_or_different_quantile[[paste0("VarBothQuantile",quantile)]] <- intersect(rownames(CtrlVarQ),rownames(HSVarQ))
  CtrlMeanQ <- subset(MeanVarTestResults, MeanCtrlQuantile == quantile)
  HSMeanQ <- subset(MeanVarTestResults, MeanHSQuantile == quantile)
  shared_or_different_quantile[[paste0("MeanCtrlQuantile",quantile)]] <- setdiff(rownames(CtrlMeanQ),rownames(HSMeanQ))
  shared_or_different_quantile[[paste0("MeanHSQuantile",quantile)]] <- setdiff(rownames(HSMeanQ),rownames(CtrlMeanQ))
  shared_or_different_quantile[[paste0("MeanBothQuantile",quantile)]] <- intersect(rownames(CtrlMeanQ),rownames(HSMeanQ))
}

View(shared_or_different_quantile)

# Run the enrichment
gProfiler_results <- gProfiler_genesets(gene_sets = shared_or_different_quantile, background_genes = background_genes)
all_results <- gProfiler_results[[1]]
top_20_results <- gProfiler_results[[2]]

# Split dataframe into list of dataframes based on the dataset
all_results_split = split(all_results, f = all_results$dataset)
top_20_results_split = split(top_20_results, f = top_20_results$dataset)

# What types of GO terms is the dataset enriched in?
unique(all_results$source)

##################################
##### (3)  GO Results ############
##################################

GO_types= c("BP","MF","CC")

for (GO_type in GO_types){

# Top 20 result plots
top20_results_plots <- sapply(names(top_20_results_split), 
                          plot_top_GO, 
                           data = top_20_results_split, 
                           GO_type = paste0("GO:",GO_type),
                           simplify = FALSE,
                           USE.NAMES = TRUE)
for (category in names(top_20_results_split)){
  ggsave(dpi = 300,plot=top20_results_plots[[category]],
        filename = paste0("GO_",GO_type,"_Top20_",category,".svg"))
}

# All results plots
all_results_GO <- sapply(names(all_results_split), 
                           subset_GO, 
                           data = all_results_split, 
                           GO_type = paste0("GO:",GO_type),
                           simplify = FALSE,
                           USE.NAMES = TRUE)

# Look at the overlaps
for (quantile in quantiles){
  GO_Venn_MeanQuant <- GO_Threeset_Venn(GO_data = all_results_GO,
                                          gene_data = shared_or_different_quantile,
                                          quantile = quantile,
                                          is_data_reduced = F,
                                          GO_category = GO_type,
                                          Var_or_Mean = "Mean")
  print(GO_Venn_MeanQuant)
  ggsave(dpi = 300,plot=GO_Venn_MeanQuant,
         filename = paste0("GO_",GO_type,"_Venn_MeanQuant",quantile,".svg"))
  
  GO_Venn_VarQuant <- GO_Threeset_Venn(GO_data = all_results_GO,
                                         gene_data = shared_or_different_quantile,
                                         quantile = quantile,
                                         is_data_reduced = F,
                                         GO_category = GO_type,
                                         Var_or_Mean = "Var")
  print(GO_Venn_VarQuant)
  ggsave(dpi = 300,plot=GO_Venn_VarQuant,
         filename = paste0("GO_",GO_type,"_Venn_VarQuant",quantile,".svg"))
}
}

##########################################
##### End) Save the Work Environment #####
##########################################

# Directory to save
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\MAD")

# Save all enrichments for the top 5% variability
GO_Table_5percentLowestVarGenes <- sanitize_for_csv(all_results_split[["VarBothQuantile1"]])
write.csv(GO_Table_5percentLowestVarGenes,file = 'GO_Table_5percentLowestVarGenes.csv',row.names = F)
GO_Table_5percentHighestVarGenes <- sanitize_for_csv(all_results_split[["VarBothQuantile20"]])
write.csv(GO_Table_5percentHighestVarGenes,file = 'GO_Table_5percentHighestVarGenes.csv',row.names = F)

# Working environemnt
save.image(file='TopBottomQuantileAnalyses_MAD.RData')
