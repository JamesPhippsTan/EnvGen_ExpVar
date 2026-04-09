Contains R, bash, Python, and D scripts used to generate figures and tables within the manuscript 'Environmental perturbation increases gene expression variability and unmasks genetic regulation for transcriptional robustness' (https://www.biorxiv.org/content/10.64898/2026.02.18.706644v1.full). The starting datasets, as well as key datasets outputted across various stages of the analysis, can be found here (https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.47SZWA).


This code repository is split into the following folders and files:

MAD - variability rank analyses using the median absolute deviation metric (MAD) (Figure 2).

DE_DV - differential expression (DE) and differential variability (DV) analysis between diets (Figure 3 and supplement).

Developmental_time - egg-adult survival eclosion day statistics (Figure 3 and supplement).

GraVe_Mapping -  analyses of genetic variation - eQTL (mean expression-regulating quantitative trait loci) and veQTL (expression variability-regulating quantitative trait loci) (Figure 1, Figure 3 and supplement). Named because it hinges on the bioinformatics tools Gra(mmar) by Yurii Aulchenko and Ve(qtl mapper) by A.A.Brown. Also contains modified versions of tensorqtl python codes (substitute these files into an existing conda installation of tensorqtl) to map eQTL and veqtl mapper D codes (requires following https://funpopgen.github.io/veqtl-mapper/ to build from source but using the files here instead of those in the default git clone) to map veQTL.

FDI_cal - used to calculate the fraction (F) of derived (D) SNPs that increase (I) either mean expression for eQTL or expression variability for veQTL and then plot distributions of this fraction across subsamples for a defined set of SNPs, e.g., those that are QTL or those that are not QTL but with the same MAF distribution as QTL (Figure 4 and supplement).

pipelines - snakemake pipelines integrating analyses and codes in GraVe_Mapping and FDI_cal made by Huiting Xu. Speeds up reproducing the results in the manuscript and facilitates adaptation to more complex experimental designs in the future.  

R scripts 0 until 6 - RNAseq data processing and quality checks (Figure 1 and supplement).

_Functions.R - functions for specific types of analyses
