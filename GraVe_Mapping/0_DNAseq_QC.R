# Check various quality metrics of the final DNAseq dataset
# Script by James Tan

# Last Updated: 10/6/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(ggplot2)
library(patchwork)
library(svglite)


#######################
##### Datasets ########
#######################

# The genome coverage per sample is obtained within the following info metadata file
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
head.info  <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)

# The samples with RNAseq AND DNAseq
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
Ctrl_flies <- read.table("Ctrl_samples.txt",h=F)
HS_flies <- read.table("HS_samples.txt",h=F)

# Site Missingness - for a SNP, which proportion of the whole population has been genotyped at this SNP
Ctrl_Site_Missingness <- read.table("Dmel_Ctrl_SiteMissingness.lmiss",header = T,sep='\t', check.names = FALSE)
HS_Site_Missingness <- read.table("Dmel_HS_SiteMissingness.lmiss",header = T,sep='\t', check.names = FALSE)

# Per-fly SNP missigness
Ctrl_Fly_Missingness  <- read.table("Dmel_Ctrl_IndivMissingness.imiss",header = T,sep='\t', check.names = FALSE)
HS_Fly_Missingness  <- read.table("Dmel_HS_IndivMissingness.imiss",header = T,sep='\t', check.names = FALSE)


#######################
##### QC Plotting #####
#######################

# Coverage over each condition
rownames(head.info) <- gsub('_','',head.info$id)

# Set up the coverage dataframe
Condition_df <- c(rep("Control", nrow(Ctrl_flies)), 
                  rep("HS", nrow(HS_flies)))
Ctrl_HS_Fly_IDs <- c(Ctrl_flies$V1,HS_flies$V1)
Coverage_df <- data.frame(Condition=Condition_df,
           GenomeCoverage=head.info[Ctrl_HS_Fly_IDs,'coverage'])
Coverage_plot <- ggplot(Coverage_df, aes(y=GenomeCoverage,x = Condition,fill=Condition)) + 
  geom_boxplot(col='darkgrey') +
  ylab(expression('Genome coverage per fly / x'))+
  xlab('Condition')+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  theme_classic()  +
  scale_fill_manual(values = c("Control" = "#C6B49E", "HS" = "#DF9F65"))
Coverage_plot

# Per-SNP genotyping rate (1-Missigness)
SNP_df <- c(rep("Control", nrow(Ctrl_Site_Missingness)), 
                  rep("HS", nrow(HS_Site_Missingness)))
GenotypingRateSNP_df <- data.frame(Condition=SNP_df,
                          GenotypingRate=1-c(Ctrl_Site_Missingness$F_MISS,
                                           HS_Site_Missingness$F_MISS))
GenotypingRateSNP_plot <- ggplot(GenotypingRateSNP_df, aes(y=GenotypingRate*100,x = Condition,fill=Condition)) + 
  geom_boxplot(col='darkgrey') +
  ylab('% of population genotyped per SNP')+
  xlab('Condition')+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  theme_classic()  +
  ylim(0,100)+
  scale_fill_manual(values = c("Control" = "#C6B49E", "HS" = "#DF9F65"))
GenotypingRateSNP_plot

# per-fly completeness
Fly_df <- c(rep("Control", nrow(Ctrl_Fly_Missingness)), 
            rep("HS", nrow(HS_Fly_Missingness)))
GenotypingRateFly_df <- data.frame(Condition=Fly_df,
                                   GenotypingRate=1-c(Ctrl_Fly_Missingness$F_MISS,
                                                      HS_Fly_Missingness$F_MISS))
GenotypingRateFly_plot <- ggplot(GenotypingRateFly_df, aes(y=GenotypingRate*100,x = Condition,fill=Condition)) + 
  geom_boxplot(col='darkgrey') +
  ylab('% of all SNPs genotyped per fly')+
  xlab('Condition')+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))+
  theme_classic()  +
  ylim(0,100)+
  scale_fill_manual(values = c("Control" = "#C6B49E", "HS" = "#DF9F65"))
GenotypingRateFly_plot


# Put all plots together
QC_plots <-  Coverage_plot+ theme(legend.position = "none") + 
  GenotypingRateSNP_plot + theme(legend.position = "none")+ 
  GenotypingRateFly_plot+ theme(legend.position = "none")+ 
  plot_layout(ncol = 3)

QC_plots

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
ggsave(plot = QC_plots,filename = "QC_plots_DNAseq.svg",height=4,width=8,dpi=300)
