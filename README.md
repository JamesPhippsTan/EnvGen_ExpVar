Contains R, python, and bash scripts utilised to generate the figures and tables within 'Environmental perturbation increases gene expression variability and unmasks genetic regulation for transcriptional robustness' (https://www.biorxiv.org/content/10.64898/2026.02.18.706644v1.full). The starting datasets, as well as intermediate datasets outputted across various stages of the analysis, can be found here (https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.47SZWA).


This code repository is split into the following fildes and folders:

MAD - variability rank analyses using the median absolute deviation metric (MAD) (Figure 2).

DE_DV - differential expression (DE) and differential variability (DV) analysis between diets (Figure 3 and supplement).

Developmental_time - egg-adult survival eclosion day statistics (Figure 3 and supplement).

GraVe_Mapping - Gra(mmar) Ve(qtl) mappping in each diet (Figure 4 and supplement). Contains modified versions of tensorqtl python codes (substitute these files into an existing conda installation) and veqtl mapper D code (requires following https://funpopgen.github.io/veqtl-mapper/ to build from source but using the files here instead of those in the default git clone). 

FDI_cal - used to calculate the fraction (F) of derived (D) SNPs that increase (I) the metric of interest - either mean expression for eQTL or variability in expression for veQTL - and then plot distributions of this fraction across subsamples for a defined set of SNPs (e.g., those that are QTL or those that are not QTL).

pipelines - snakemake pipelines integrating analyses and codes in GraVe_Mapping and FDI_cal. Speed up reproducing the results in the manuscript and facilitates adaptation to more complex experimental designs.  

R scripts 0) until 7) - RNAseq data processing and quality checks (Figure 1 and supplement).

_Functions.R - functions for specific but inter-related aspects of the analyses
