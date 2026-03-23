Contains R, python, and bash scripts utilised to generate the figures and tables within 'Environmental perturbation increases gene expression variability and unmasks genetic regulation for transcriptional robustness' (https://www.biorxiv.org/content/10.64898/2026.02.18.706644v1.full). It is split into the following subfolders:


Start page - helper functions and RNAseq data QC (Figure 1 and supplement).

MAD - contains variability rank analyses using the median absolute deviation metric (MAD) (Figure 2).

DE_DV - differential expression (DE) and differential variability (DV) analysis between diets (Figure 3 and supplement).

Developmental_time - used for egg-adult survival eclosion day statistics (Figure 3 and supplement).

GraVe_Mapping - Gra(mmar) Ve(qtl) mappping in each diet (Figure 4 and supplement). Contains modified versions of tensorqtl python codes (substitute these files into an existing conda installation) and veqtl mapper D code (requires following https://funpopgen.github.io/veqtl-mapper/ to build from source but using the files here instead of those in the default git clone). Snakemake pipeline in progress for better flow and ease of reproducibility.

FDI_cal - used to calculate the fraction (F) of derived (D) SNPs that increase (I) the metric of interest, either mean expression for eQTL or variability in expression for veQTL, and plot a distribution of this fraction for a defined set of SNPs. 

pipelines - snakemake pipelines integrating analyses and codes in GraVe_Mapping and FDI_cal to allow (1) quicker reproducibility of the results in the manuscript and (2) adaptation to more complex experimental designs.  
