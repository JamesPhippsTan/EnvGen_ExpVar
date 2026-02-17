# PCA of the genotypes in the veQTL mapping population
# Last updated: 4/5/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

## Load 
library('snpStats')
library(vcfR)
library("gdsfmt")
library("SNPRelate")
library(ggplot2)

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
load(file='0a_PopulationPCA.RData')

#######################
##### Datasets ########
#######################

# Read conditions
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
cov <- read.table('Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt',header = T)
environment <- data.frame(sample.id=gsub("_", "", cov$id),environment=as.factor(cov$treatment))

# Read the fly genotype matrix into R
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
bed.fn <- './Dmel_Ctrl_HS_final_full.bed'
fam.fn <- './Dmel_Ctrl_HS_final_full.fam'
bim.fn <- './Dmel_Ctrl_HS_final_full.bim'

# Create merged genofile
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "Dmel_Ctrl_HS_final_full.gds",cvt.chr="char")

snpgdsSummary("Dmel_Ctrl_HS_final_full.gds")


###########################
##### PCA plotting ########
###########################

genofile <- snpgdsOpen("Dmel_Ctrl_HS_final_full.gds")

# Make PCs with all 413348 SNPs (they have been LD-pruned already)
pca <- snpgdsPCA(genofile, autosome.only=F, num.thread=4)
pc.percent <- pca$varprop*100
round(pc.percent, 2)[1]

# Plot the first two PCSs by environment
PC1and2 <- data.frame(sample.id = gsub("-", "_",pca$sample.id),
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
PC1and2_env <- merge(PC1and2,environment,by='sample.id')
PC1and2_env$Diet <- ifelse(PC1and2_env$environment==1,'Control','High Sugar')


# Plot
PCplot <- ggplot(data=PC1and2_env,aes(x=EV1,y=EV2,group=Diet,col=Diet))+geom_point()+ 
  scale_color_manual(values=c('#C6B49F','#DF9F65'))+ 
  ylab(paste0('SNPs PC2 (',round(pc.percent, 2)[2],'%)')) + 
  xlab(paste0('SNPs PC1 (',round(pc.percent, 2)[1],'%)')) +
  theme_classic()+  # Clean theme with larger font

theme(
    legend.position = c(0.15, 0.85),  # Places legend inside plot (x,y between 0-1)
    legend.background = element_rect(fill = "white", color = "black", size = 0.1),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) 

PCplot_nolegend <- PCplot+ theme(legend.position = "none")
PCplot_nolegend
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
ggsave(filename = 'PopulationPCA.png',plot = PCplot,height = 3,width = 3)
ggsave(filename = 'PopulationPCAnolegend.png',plot = PCplot_nolegend,height = 3,width = 3)

#
# Save working environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping")
save.image(file='0a_PopulationPCA.RData')
