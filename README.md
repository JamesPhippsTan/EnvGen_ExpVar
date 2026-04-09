Contains R, python, and bash scripts utilised to generate the figures and tables within 'Environmental perturbation increases gene expression variability and unmasks genetic regulation for transcriptional robustness' (https://www.biorxiv.org/content/10.64898/2026.02.18.706644v1.full). The starting datasets, as well as key datasets outputted across various stages of the analysis, can be found here (https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.47SZWA).


This code repository is split into the following folders and files:

MAD - variability rank analyses using the median absolute deviation metric (MAD) (Figure 2).

DE_DV - differential expression (DE) and differential variability (DV) analysis between diets (Figure 3 and supplement).

Developmental_time - egg-adult survival eclosion day statistics (Figure 3 and supplement).

GraVe_Mapping - Gra(mmar) Ve(qtl) mappping in each diet (Figure 4 and supplement). Also contains modified versions of tensorqtl python codes (substitute these files into an existing conda installation of tensorqtl) to map eQTL - mean expression-regulating genetic loci - and veqtl mapper D codes (requires following https://funpopgen.github.io/veqtl-mapper/ to build from source but using the files here instead of those in the default git clone) to map veQTL - expression variability-regulating genetic loci.

FDI_cal - used to calculate the fraction (F) of derived (D) SNPs that increase (I) either mean expression for eQTL or expression variability for veQTL and then plot distributions of this fraction across subsamples for a defined set of SNPs (e.g., those that are QTL or those that are not QTL but with the same MAF distribution as QTL).

pipelines - snakemake pipelines integrating analyses and codes in GraVe_Mapping and FDI_cal made by Huiting Xu. Speeds up reproducing the results in the manuscript and facilitates adaptation to more complex experimental designs in the future.  

R scripts 0 until 6 - RNAseq data processing and quality checks (Figure 1 and supplement).

_Functions.R - functions for specific types of analyses
