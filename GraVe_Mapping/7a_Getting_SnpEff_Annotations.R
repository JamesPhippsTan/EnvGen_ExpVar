# Getting SnpEff annotations
# For use with other SNP regulatory feature enrichments

# Last Updated: 7/2/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls(all = T))

library(tidyr)
library(dplyr)

# Load saved script environment
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
load(file='7a_Getting_SnpEff_Annotations')

#######################
##### Datasets ########
#######################

# Position file to re-map actual SNP names to trans-veQTL
setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\")
SnpEff_Annotations_Full <- read.table("SnpEff_Annotated_SNP_Set.txt",header = F)

# Condense to ID and Annotations
SnpEff_Annotations <- SnpEff_Annotations_Full[,c("V3","V8")]
colnames(SnpEff_Annotations) <- c('SNP',"Annotation")
View(SnpEff_Annotations)

# Create a list of features I want
SnpEff_features <- c("exon_variant",
                     "start_lost",
                     "stop_gained",
                     "stop_lost",
                     "missense_variant",
                     "synonymous_variant",
                     "3_prime_UTR_variant",
                     "5_prime_UTR_variant",
                     "intron_variant",
                     "intergenic_region",
                     "splice_region_variant")

# For each feature, find them within the annotation column
# Then, make a column indicating whether the variant has that annotation
for (feature in SnpEff_features) {
  col_name <- feature
  SnpEff_Annotations <- SnpEff_Annotations %>%
    mutate(!!col_name := ifelse(str_detect(Annotation, feature), "yes", "no"))
}

colnames(SnpEff_Annotations) <- gsub("_variant","",colnames(SnpEff_Annotations))

# Remove the loooong annotations
SnpEff_Annotations_final <- SnpEff_Annotations[,-2]

###################################################
##### (End) Save working environment ##############
###################################################

write.csv(SnpEff_Annotations_final,"SnpEff_Annotations_final.csv",row.names = F)

save.image(file='7a_Getting_SnpEff_Annotations')
