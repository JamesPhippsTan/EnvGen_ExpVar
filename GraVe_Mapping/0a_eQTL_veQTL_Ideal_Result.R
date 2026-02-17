# Code to plot the difference between a hypothetical eQTL and a veQTL
# Script by James Tan

# Last updated: 21/5/2025


#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(ggplot2)

# Working Directory 
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\Ideal")

#######################
##### Datasets ########
#######################

# Metadata on all samples, e.g., conditions and batches
Ideal_dataset   <- read.csv("Ideal_eQTL_veQTL_results.csv",h=T)

# Neither
Neither <- subset(Ideal_dataset, Include_in_which_plot == 'Neither')
Neitherplot <- ggplot(Neither, aes(x=Genotype.at.Genomic.Site, y=Expression)) + 
  geom_boxplot(color="darkgrey", fill="white", alpha=0.2) + 
  labs(x = "Genotype") + 
  ylim(5,15)+ 
  labs(y = "Transcript Level") + theme_classic() + theme(axis.text.y = element_blank()) + theme(axis.ticks.y = element_blank())
Neitherplot

# eQTLs
eQTL <- subset(Ideal_dataset, Include_in_which_plot == 'eQTL')
eQTLplot <- ggplot(eQTL, aes(x=Genotype.at.Genomic.Site, y=Expression)) + 
  geom_boxplot(color="darkgrey", fill='#4F8E4D', alpha=0.2) + 
  labs(x = "eQTL genotype") + 
  labs(y = "Transcript Level") + theme_classic() + theme(axis.text.y = element_blank()) + theme(axis.ticks.y = element_blank())
eQTLplot

# veQTLs
veQTL <- subset(Ideal_dataset, Include_in_which_plot == 'veQTL')
veQTLplot <- ggplot(veQTL, aes(x=Genotype.at.Genomic.Site, y=Expression)) + 
  geom_boxplot(color="darkgrey", fill="#611BB8", alpha=0.2) + 
  labs(x = "veQTL genotype") + 
  labs(y = "Transcript Level") + theme_classic() + theme(axis.text.y = element_blank()) + theme(axis.ticks.y = element_blank())
veQTLplot

# Both
Both <- subset(Ideal_dataset, Include_in_which_plot == 'Both')
Bothplot <- ggplot(Both, aes(x=Genotype.at.Genomic.Site, y=Expression)) + 
  geom_boxplot(color="darkgrey", fill="red", alpha=0.2) + 
  labs(x = "Genotype") + 
  labs(y = "Transcript Level") + theme_classic() + theme(axis.text.y = element_blank()) + theme(axis.ticks.y = element_blank())
Bothplot

#######################
##### Save plots ######
#######################

png(file = 'example_eQTL.png',  width = 2, height = 2,units = 'in',res = 144)
eQTLplot
dev.off()

png(file = 'example_veQTL.png',  width = 2, height = 2,units = 'in',res = 144)
veQTLplot
dev.off()

png(file = 'example_neithereQTLveQTL.png',  width = 2, height = 2,units = 'in',res = 144)
Neitherplot
dev.off()

png(file = 'example_botheQTLveQTL.png',  width = 2, height = 2,units = 'in',res = 144)
Bothplot
dev.off()
