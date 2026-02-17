# Applying the GRAMMAR correction - removing GRM from the data - Ctrl dataset

# Last updated: 7/11/2024

#################################
##### Packages and Setup ########
#################################

rm(list = ls())
library('GenABEL')

setwd("/tmp/global2/jtanshengyi/veQTL_GraVe/data")

# Get the number of the specific array job
array_job_number <- as.integer(Sys.getenv("SGE_TASK_ID"))

# Define the number of genes per job 
# As we have 8763 genes and I have setup 40 jobs, the number of genes is around 220
n_genes_per_job=220

##################################################
##### GRAMMAR to remove covariates and GRM  ######
##################################################

# Take the Genotype and Phenotype files and convert to a format operable by GenAbel
GP_Data <- load.gwaa.data(genofile = 'Dmel_Ctrl_MAF5_Miss50_GRAMMAR.raw',phenofile = 'GraVe_phenotypes_Ctrl.txt',id='id')

# GRM to estimate unknown relatedness - all individuals, all SNPs
GRM <- ibs(GP_Data, w="freq")

# Retrieve the subset of phenotypes to be transformed 
# Begin with
start <- (array_job_number-1)*n_genes_per_job+3
end <- start+n_genes_per_job-1

# It is possible that the 'end' index may go beyond the last column and thus not exist. In this case, recode 'end' as the index of the last column
if (end > ncol(GP_Data@phdata)) {
  end <- ncol(GP_Data@phdata)
}
P_Data_GRAMMARed <- GP_Data@phdata[,c(1,2,start:end)]

# For each phenotype in the subset, apply the polygenic() function to obtain GRAMMAR-corrected counts (i.e., residuals)
position = 3
for (phenotype_index in start:end) {
  polyLMM <- polygenic(GP_Data@phdata[,phenotype_index],GRM, GP_Data)
  P_Data_GRAMMARed[,position] <- polyLMM$grresidualY
  position = position + 1
}


# Get rid of the 'sex' column (it is unneeded for steps after this).
P_Data_GRAMMARed <- P_Data_GRAMMARed[,-2]

####################################
##### Save the final matrix ########
####################################

setwd("GRAMMARed_phenotypes_Ctrl")
first_gene <- start-2
last_gene <- end-2
write.table(P_Data_GRAMMARed,paste0('GRAMMARed_phenotypes_Ctrl_',first_gene,'_',last_gene,'.txt'),quote=F)


