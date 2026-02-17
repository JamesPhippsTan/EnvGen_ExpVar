# Concatenating the various downsamplings
# Running them together takes way too long
# Run them in parallel then merge

# Last Updated: 11/9/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(DESeq2)
library(dplyr)
library(edgeR)
library(limma)
library(tidyr)
library(tibble)
library(car)
library(matrixStats)
library(ggplot2)
library(ggrepel)
library(gamlss)

####################
##### Datasets #####
####################

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_Downsampled")
load(file='DE_DV_Downsampled_Series_50100.RData')
downsampling_DE_df_full <- downsampling_DE_df 
downsampling_DV_df_full <- downsampling_DV_df

subsampling_results <- c('DE_DV_Downsampled_Series_Small_200300.RData',
                      'DE_DV_Downsampled_Series_Small_500.RData',
                      'DE_DV_Downsampled_Series_Small_900.RData')

for (result in subsampling_results){
  load(file=result)
  downsampling_DE_df_full <- rbind(downsampling_DE_df_full,downsampling_DE_df) 
  downsampling_DV_df_full <- rbind(downsampling_DV_df_full,downsampling_DV_df)
}

# Get the correct levels
correct_levels <- c("50", "100", "200", "300", "500", "900", "Full")
downsampling_DE_df_full$n <- factor(downsampling_DE_df_full$n, levels = correct_levels)
downsampling_DV_df_full$n <- factor(downsampling_DV_df_full$n, levels = correct_levels)

############################
##### Plot the results #####
############################

# Plot number of DEGs
DEG_number_by_subsample_size <- ggplot(data = downsampling_DE_df_full, aes(x=n,y=DEGs)) + 
  geom_boxplot(fill='#4F8E4D',col='darkgrey')+
  geom_point(col='darkgrey')+
  theme_classic()+  
  ylab('Number of DEGs')+
  xlab('Sample size per condition')+  
  ylim(0,8763)+
  geom_hline(yintercept=8763, col="red")
DEG_number_by_subsample_size

DVG_number_by_subsample_size <- ggplot(data = downsampling_DV_df_full, aes(x=n,y=DVGs)) + 
  geom_boxplot(fill='#611BB8',col='darkgrey')+
  geom_point(col='darkgrey')+
  theme_classic()+
  ylab('Number of DVGs')+
  xlab('Sample size per condition')+
  ylim(0,8763)+
  geom_hline(yintercept=8763, col="red")
DVG_number_by_subsample_size

# As a fraction of the full sample size
Fraction_DEG_by_subsample_size <- ggplot(data = downsampling_DE_df_full, aes(x=n,y=100*DEGs/8763)) + 
  geom_boxplot(fill='#4F8E4D',col='darkgrey')+
  geom_point(col='darkgrey')+
  theme_classic()+  
  ylab('% of transcriptome DE')+
  xlab('Sample size per condition')+  
  ylim(0,100)
Fraction_DEG_by_subsample_size

Fraction_DVG_by_subsample_size <- ggplot(data = downsampling_DV_df_full, aes(x=n,y=100*DVGs/8763)) + 
  geom_boxplot(fill='#611BB8',col='darkgrey')+
  geom_point(col='darkgrey')+
  theme_classic()+
  ylab('% of transcriptome DV')+
  xlab('Sample size per condition')+
  ylim(0,100)
Fraction_DVG_by_subsample_size

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_Downsampled")
ggsave(dpi = 300,plot = Fraction_DEG_by_subsample_size,filename = "DEG_fraction_by_subsample_size.svg",width = 3,height = 3)
ggsave(dpi = 300,plot = Fraction_DVG_by_subsample_size,filename = "DVG_fraction_by_subsample_size.svg",width = 3,height = 3)

ggsave(dpi = 300,plot = DEG_number_by_subsample_size,filename = "DEG_number_by_subsample_size.svg",width = 3,height = 3)
ggsave(dpi = 300,plot = DVG_number_by_subsample_size,filename = "DVG_number_by_subsample_size.svg",width = 3,height = 3)



##################################
##### End) Save the workspace ####
##################################

# Save workspace
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\DE_DV\\DE_DV_Downsampled/")
save.image(file='DE_DV_Downsampled_Series_Concatenated.RData')
