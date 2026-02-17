# Turn the ped and map files into raw files readable by GenAbel

# Last updated: 7/11/2024

#################################
##### Packages and Setup ########
#################################

rm(list = ls())
library('GenABEL')

setwd("/tmp/global2/jtanshengyi/veQTL_GraVe/data")

convert.snp.ped('Dmel_Ctrl_MAF5_Miss50_GRAMMAR.ped','Dmel_Ctrl_MAF5_Miss50_GRAMMAR.map','Dmel_Ctrl_MAF5_Miss50_GRAMMAR.raw')
convert.snp.ped('Dmel_HS_MAF5_Miss50_GRAMMAR.ped','Dmel_HS_MAF5_Miss50_GRAMMAR.map','Dmel_HS_MAF5_Miss50_GRAMMAR.raw')



