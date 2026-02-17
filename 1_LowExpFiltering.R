# Filtering out low expression genes and quality checking

# Original title:
# DEC 29th 2020
# GETTING GENE EXPRESSION MATRIX FROM HEAD that includes both conditions HS and CTRL
# for GxE mapping

# Last Updated: 6/6/25

#################################
##### Packages and Setup ########
#################################


rm(list = ls(all = T))
library(dplyr)

######################
##### Dataset ########
######################

#Head
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")
head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)
head.info$plate <- gsub("b", "", head.info$plate)
#head.counts <- read.table("UnprocessedCountMatrix.txt",h=T,check.names = F)

dim(head.counts) #17727 genes, 2136 samples

# FIRST: FILTER OUT LOW COUNT GENES BASED ON MEAN CPM <1 
{
  #Transform to CPM to adjust for library size before doing any filtering, do this with both diets together
  
  #Head
  library(edgeR)
  head.cpm <- DGEList(counts = head.counts)
  head.cpm <- calcNormFactors(object = head.cpm)
  head.cpm <- cpm(head.cpm, log = F)  
  
  # Look at the library sizes
  library_sizes <- head.cpm$samples$lib.size*head.cpm$samples$norm.factors
  
  
  # Get control samples (n=989)  
  head.c <- head.counts[,colnames(head.counts) %in% head.info$id[head.info$treatment==1]]
  head.cpm.c <- head.cpm[ , colnames(head.cpm) %in% colnames(head.c)] ; dim(head.cpm.c)
  
  # Get hs samples (n=1147)  
  head.hs <- head.counts[,colnames(head.counts) %in% head.info$id[head.info$treatment==2]]
  head.cpm.hs <- head.cpm[ , colnames(head.cpm) %in% colnames(head.hs)]; dim(head.cpm.hs)
  
  #keep genes with average CPM > 1 and CPM>1 in at least 20% of the individuals, per treatment
  #(keeps = 8800)
  phenotyping.rate.c <- apply(head.cpm.c, 1, function(x) sum(x>1)/length(x))
  head.cpm1.c <- subset(head.cpm.c, (rowSums(head.cpm.c)/ncol(head.cpm.c))>1 & phenotyping.rate.c >0.2) ; dim(head.cpm.c) ; dim(head.cpm1.c)
  phenotyping.rate1.c <- apply(head.cpm1.c, 1, function(x) sum(x>1)/length(x))
  
  #(keeps = 9707) #quite a diff with control
  phenotyping.rate.hs <- apply(head.cpm.hs, 1, function(x) sum(x>1)/length(x))
  head.cpm1.hs <- subset(head.cpm.hs, (rowSums(head.cpm.hs)/ncol(head.cpm.hs))>1 & phenotyping.rate.hs >0.2) ; dim(head.cpm.hs) ; dim(head.cpm1.hs)
  phenotyping.rate1.hs <- apply(head.cpm1.hs, 1, function(x) sum(x>1)/length(x))
  
  #plots 
  phenotyping.rate=phenotyping.rate.hs
  phenotyping.rate1=phenotyping.rate1.hs
  head.cpm1.tmp=head.cpm1.hs
  head.cpm.tmp=head.cpm.hs
  head.counts.tmp=head.hs #put this back
  
  par(mfrow=c(2,2))
  hist(phenotyping.rate, main="%individuals CPM>1 per SNP", col="grey"); hist(phenotyping.rate1, col=adjustcolor("red",0.5), add=T)
  samplesize.site <- apply(head.cpm1.tmp, 1, function(x) sum(x>1)) #min # ind per snp in ctrl 198, en hs also 230
  hist(samplesize.site, co="grey", main="Sample size per site") 
  plot(log2(rowSums(head.cpm.tmp)/ncol(head.cpm.tmp)) ~ phenotyping.rate, ylab="average CPM per site (log2)", xlab="phenotyping rate")  ; abline(h=0, lty=2)
  points(log2(rowSums(head.cpm1.tmp)/ncol(head.cpm1.tmp)) ~ phenotyping.rate1, col="salmon")
  c <- rowSums(head.cpm.tmp)/ncol(head.cpm.tmp)
  co <- rowSums(head.counts.tmp)/ncol(head.counts.tmp)
  plot (log2(c) ~ log2(co), xlab ="Average gene counts (log2)", ylab="Average CPM per site (log2)"); abline(h=0, lty=2) ; abline(v=2, lty=2)
  
  hist(log2(c), main="average CPM per gene before filtering", col="grey", breaks=50) 
  hist(log2(rowSums(head.cpm1.tmp)/ncol(head.cpm1.tmp)), breaks=50,col=adjustcolor("salmon",0.4),add=T, main="average CPM per gene after removing low count genes")
  missdata.indiv <- apply(head.cpm.tmp, 2, function(x) sum(x>1)/length(x))
  missdata.indiv1 <- apply(head.cpm1.tmp, 2, function(x) sum(x>1)/length(x))
  hist(missdata.indiv, col="grey", main="% sites CPM >1 per individual before filter"); hist(missdata.indiv1, col="grey", add=F, main="% sites CPM >1 per individual after filter")
}      

# SECOND: OVERLAP BETWEEN GENES EXPRESSED IN THE hs AND IN ctrl
{      
  # 8786 genes overlap between the two datasets
  # 14 genes are only expressed in ctrl
  # 921 genes are only expressed in hs
  both <- intersect(rownames(head.cpm1.c), rownames(head.cpm1.hs)) ; length(both)
  length(setdiff(rownames(head.cpm1.c), rownames(head.cpm1.hs)))
  length(setdiff(rownames(head.cpm1.hs), rownames(head.cpm1.c)))
} 

# WRITE THE FROZEN VERSION OF THE FILTERED RAW GENE COUNT MATRIX - USE THIS FOR ANY FURTHER ANALYSIS
# 8786 genes expressed in both diets after filters above 
{  
  #Raw counts, but genes with CPM < 1 have been removed
  head.counts <- subset(head.counts, rownames(head.counts) %in% both)
  all.equal(colnames(head.counts), as.character(head.info$id))
  write.table(head.counts, "RawCounts_CPM1_head_hsctrl_Dec29.20.txt", col.names = T, row.names = T, sep="\t", quote=F)
    write.table(head.info, "Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt", col.names = T, row.names = F, sep="\t", quote=F)
}
