# DEC 29th 2020
# GETTING GENE EXPRESSION MATRIX FROM HEAD that includes both conditions HS and CTRL
# for GxE mapping

library(dplyr)

setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data")

# LOAD RAW GENE EXPRESSION MATRIX FOR HEAD 

  #Head
  head.info   <- read.table("Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt",h=T)
  head.info$plate <- gsub("b", "", head.info$plate)
  head.counts <- read.table("RawCounts_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt",h=T,check.names = F)

  dim(head.counts) #17727 genes, 2136 samples


# FIRST: FILTER OUT LOW COUNT GENES BASED ON MEAN CPM <1 
{
  #Transform to CPM to adjust for library size before doing any filtering, do this with both diets together
  
  #Head
  library(edgeR)
  head.cpm <- DGEList(counts = head.counts)
  head.cpm <- calcNormFactors(object = head.cpm)
  head.cpm <- cpm(head.cpm, log = F)  
  rownames(head.cpm) <- rownames(head.counts)
  
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
  
  
  # --> to do, GO terms for the genes uniquely expressed in head or body
  
  
  
  # TO DO: Relationship between genes expressed in both diets
  {
  
  mean.expression.head <- rowSums(head.cpm)/ncol(head.cpm)
  mean.expression.body <- rowSums(body.cpm)/ncol(body.cpm)
  mean.expression.head1 <- rowSums(head.cpm1[rownames(head.cpm1) %in% both,])/ncol(head.cpm1)
  mean.expression.body1 <- rowSums(body.cpm1[rownames(body.cpm1) %in% both,])/ncol(body.cpm1)
  
  par(mfrow=c(1,2))
  plot(log2(mean.expression.head) ~ log2(mean.expression.body), xlab="mean CPM per gene (log2)- BODY", ylab =" mean CPM per gene (log2) - HEAD", col="lightgrey")
  points(log2(mean.expression.head[rownames(head.cpm1)]) ~ log2(mean.expression.body[rownames(head.cpm1)]), col="red")
  points(log2(mean.expression.head[rownames(body.cpm1)]) ~ log2(mean.expression.body[rownames(body.cpm1)]), col="salmon")
  points(log2(mean.expression.head1) ~ log2(mean.expression.body1), col="darkorange")
  abline(0,1, lty=2)
  points(log2(mean.expression.head[names(mean.expression.head)=="FBgn0038484"]) ~ log2(mean.expression.body[names(mean.expression.body)=="FBgn0038484"]), col="black", pch=8)
  
  #correlation of the whole transcriptome
  text(-9,15, "rho = 0.9, r2 = 0.96")
  text(-9,13, "p = 2e-16")
  cor.test(mean.expression.head , mean.expression.body, method = "spearman")
  cor.test(mean.expression.head , mean.expression.body, method = "pearson")
  
  
  #correlation between the genes that are expressed (CPM filter) and are found in both tissues
  #Rank correlation 0.53, pearson = 0.91, p=2e-16
  plot(log2(mean.expression.head1) ~ log2(mean.expression.body1), col="darkorange", xlab="Average CPM (log2) - BODY", ylab="Average CPM (log2) - HEAD")
  abline(0,1, lty=2)
  text(5,15, "rho = 0.51, r2 = 0.97")
  text(5,14, "p = 2e-16")
  cor.test(mean.expression.head1 , mean.expression.body1, method = "spearman")
  cor.test(mean.expression.head1 , mean.expression.body1, method = "pearson")
  
}      

} 
  
# WRITE THE FROZEN VERSION OF THE FILTERED RAW GENE COUNT MATRIX - USE THIS FOR ANY FURTHER ANALYSIS
# 8786 genes expressed in both diets after filters above 
{  
 #Raw counts, but genes with CPM < 1 have been removed
head.counts <- subset(head.counts, rownames(head.counts) %in% both)
all.equal(colnames(head.counts), as.character(head.info$id))
write.table(head.counts, "C:/Users/lfpal/Documents/DecanalizationProject/GeneExpressionMatrices/head/RawCounts_CPM1_head_hsctrl_Dec29.20.txt", col.names = T, row.names = T, sep="\t", quote=F)

head.cpm1 <- subset(head.cpm, rownames(head.cpm) %in% both )
all.equal(colnames(head.cpm1), as.character(head.info$id))
write.table(head.cpm1, "C:/Users/lfpal/Documents/DecanalizationProject/GeneExpressionMatrices/head/CPMcounts_CPM1_head_hsctrl_Dec29.20.txt", col.names = T, row.names = T, sep="\t", quote=F)

write.table(head.info, "C:/Users/lfpal/Documents/DecanalizationProject/GeneExpressionMatrices/head/Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt", col.names = T, row.names = F, sep="\t", quote=F)
}

# TRANSFORM DATA FOR EQTL MAPPING USING VARIANCE STABILIZING TRANSFORMATION FROM VOOM LIMMA
{
# generates log2CPM adjusted by TMM. Data  
# to get log2CPM it adds 0.5 to each count, so there are never zeros in the final matrix
# and, because log2CPM is divided by library size (TMM), samples with zeros, but differen tlibrary sizes
# will have a negative number(log2 of <1), but it wont be the same. in VST, zeros are fixed to a number
# that's why the data doesn't look normal at the end when there are too many zeros
library(limma)

head.list <- DGEList(counts=head.counts)
head.norm <- calcNormFactors(head.list) #gets norm. factors based on TMM (controlling not for lib size, but also composition)
head.voom <- voom(head.norm, design = NULL, plot=T) #the E matrix of normalized counts is the same as in voomWithQualityWeights
# log2(counts+0.5/normalized lib.size in Millions), norm library size = totalreads*norm.factors
write.table(head.voom$E, "C:/Users/lfpal/Documents/DecanalizationProject/GeneExpressionMatrices/head/VOOMCounts_CPM1_head_hsctrl_Dec29.20.txt", col.names = T, row.names = T, 
            quote=F, sep="\t")
}


# GENERATE A COVARIATE-FREE MATRIX - 
# if this works, and covfree mapping is similar to full model, then re do this by estimating svs from the VOOM dataset
# and using the right genes
{

### jan 5 2020
### generating a cov free matrix to compare to the full model in GEMMa, because running a gxe model w/o covs runs in ~3h
### while running the full fmodel with 7 covariates takes ~12h. 
### the covs were estimated on VST transformed data, but in geema im using voom counts
### so here, i should create voom-covfee-dataset to match what's going on in GEMMA
library(limma)

counts <- read.table("C:/Users/lfpal/Documents/DecanalizationProject/GeneExpressionMatrices/head/VOOMCounts_CPM1_head_hsctrl_Dec29.20.txt", h=T, check.names = F)

cov <- read.table("Covariates_forMatrixeQTL_Mar17.txt",h=T, check.names = F)
cov2 <- t(cov)
cov2 <- as.data.frame(cov2)
cov2 <- cov2[order(match(rownames(cov2), colnames(counts))),]
all.equal(rownames(cov2), colnames(counts))
cov2$treatment <- as.factor(cov2$treatment)
cov2$RNAlibBatch <- as.factor(cov2$RNAlibBatch)
cov2$RNAseqBatch <- as.factor(cov2$RNAseqBatch)
cov2$egglayBatch<- as.factor(cov2$egglayBatch)
cov2$platingBatch <- as.factor(cov2$platingBatch)



nullmodel<-model.matrix(~cov2$treatment)

# Option 1: Adjusting counts by 1 batch effect at a time (i.e., batch effects not modelled as covariates).

#adjust for batch effects that I'm using in the full GEMMA model, can only fit two batches at a time, and multiple covs. 
#when multiple batches and covs, all considered as additive effects.
#the order in which I account for batches matters, so trying to fit them all in a covariate file
#the order matters
counts.new <- limma::removeBatchEffect(counts, batch=cov2$RNAlibBatch, design = nullmodel) 
counts.new1.2 <- limma::removeBatchEffect(counts.new, batch=cov2$egglayBatch, design = nullmodel) 

counts.new2 <- limma::removeBatchEffect(counts, batch=cov2$egglayBatch, design = nullmodel) 
counts.new2.2 <- limma::removeBatchEffect(counts.new2, batch=cov2$RNAlibBatch, design = nullmodel) 

# Option 2: Adjusting counts by many batch effects simultaneously using the covariate option of removeBatchEffects

# Test 1 - just one batch effect
#one cov file including all batches factors as numerics in the covariate 
#if i just fit it using the normal contrasts set up in R which is contr.treatment
#I get different results from running it as a factor
#I guess because anovalike functions in R like a lm in limma uses contr.sum instead of contr.treatment
#setting the base line as the last level, and confining the levels to add up to 1
covariates <- model.matrix(~cov2$RNAlibBatch) 
counts.new3 <- limma::removeBatchEffect(counts, covariates = as.numeric(cov2$RNAlibBatch) ,  design = nullmodel) 

# Test 2 - just one batch effect but changing type of contrast matrix to contrast.sum() rather than the default contrast.treatment()
#oone cov file including all batches factors as numerics in the covariate 
#if i set up the contrasts to contr.sum, then I get the same results using the batch as factor or numeric,
#and then I can fit all batches as covaraites which will be similar to the actual lm in gemma
#im resetting the contrasts on a new variable not to mess up with the actual structure of the data
#following: https://support.bioconductor.org/p/85202/
#http://faculty.nps.edu/sebuttre/home/r/contrasts.html
rnalib <- cov2$RNAlibBatch
contrasts(rnalib) <- contr.sum(levels(rnalib))
covariates <- model.matrix(~rnalib)
counts.new4 <- limma::removeBatchEffect(counts, covariates = covariates[,-1],  design = nullmodel) 

## ok, now running it with all batches as covariates
rnalib <- cov2$RNAlibBatch; rnaseq <- cov2$RNAseqBatch; egg<-cov2$egglayBatch; plate<-cov2$platingBatch
contrasts(rnalib) <- contr.sum(levels(rnalib)); contrasts(rnaseq) <- contr.sum(levels(rnaseq));contrasts(egg) <- contr.sum(levels(egg)); contrasts(plate) <- contr.sum(levels(plate))
covariates <- model.matrix(~rnalib+rnaseq+egg+plate+cov2$sv1+cov2$sv2)
counts.new5 <- limma::removeBatchEffect(counts, covariates = covariates[,-1],  design = nullmodel) 

write.table(counts.new5, "C:/Users/lfpal/Documents/DecanalizationProject/GeneExpressionMatrices/head/VOOMCounts_CPM1_head_hsctrl_covfree_Jan5.21.txt",col.names = T,row.names = T,sep="\t",quote=F)

}

