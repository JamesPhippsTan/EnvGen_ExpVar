#!/usr/bin/env Rscript

# Author: Huiting
library(ggplot2)
library(svglite)

# calculate ecdf
ecdf_calculation <- function(df, file_name, true_value) {
  F <- ecdf(df$FDD)# Compute percentile
  percentile <- F(true_value)
  final_percentile <- 1-percentile
  print(file_name)
  print(final_percentile)
}

# plot qtl
qtl_plot <- function(file_path, name, true_value,color) {
  df <- read.csv(file_path)
  ecdf_calculation(df,name, true_value)
  label_hjust <- ifelse(true_value < mean(c(0.1, 0.8)), -0.05, 1.05)#(0.45,0.65)
  
  ggplot(df, aes(x = FDD)) +
    geom_histogram(
      aes(y = after_stat((count / sum(count))*100)),
      breaks = seq(0.1,0.8, by = 0.005),
      position = "identity",
      color = "black",
      fill = color
    ) +
    geom_vline(
      xintercept = true_value,
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    #annotate("text", x = true_value, y = 35, label = name, hjust = label_hjust, size = 2) +
    coord_cartesian(xlim = c(0.1, 0.8), ylim = c(0, 35)) +
    labs(
      x = "Fraction minor alleles with decreasing effect",
      y = "Percent"
    ) +
    theme_classic()
}

# Minor
# trans FMD average corrected on 24/09/25
p1 <- qtl_plot("veqtl_trans_v2/Minor_map_nonveqtl/Ctrl_FDD_summary.csv", "Ctrl trans-veQTL", 0.73,"#611BB8") # threshold is the FDD without between
p2 <- qtl_plot("veqtl_trans_v2/Minor_map_nonveqtl/HS_FDD_summary.csv", "HS trans-veQTL", 0.8970998,"#611BB8") # threshold is the average FDD of veqtl SNPs
p3 <- qtl_plot("eqtl_trans/Minor_map_noneqtl/Ctrl_FDD_summary.csv", "Ctrl trans-eQTL", 0.599386,"#4F8E4D")# threshold is the average FDD of eqtl SNPs
p4 <- qtl_plot("eqtl_trans/Minor_map_noneqtl/HS_FDD_summary.csv", "HS trans-eQTL", 0.5811372,"#4F8E4D")# threshold is the average FDD of eqtl SNPs

ggsave("Minor_veQTL_Ctrl_FDD_plot.svg", p1, width = 4, height = 4,dpi=300)
ggsave("Minor_veQTL_HS_FDD_plot.svg", p2, width = 4, height = 4,dpi=300)
ggsave("Minor_eQTL_Ctrl_FDD_plot.svg", p3, width = 4, height = 4,dpi=300)
ggsave("Minor_eQTL_HS_FDD_plot.svg", p4, width = 4, height = 4,dpi=300)

"""
# trans
p1 <- qtl_plot("veqtl_trans_v2/MAF_map_nonveqtl/Ctrl_FDD_summary.csv", "Ctrl trans-veQTL", 0.5776124853843796,"#611BB8") # threshold is the FDD without between
p2 <- qtl_plot("veqtl_trans_v2/MAF_map_nonveqtl/HS_FDD_summary.csv", "HS trans-veQTL", 0.6558296,"#611BB8") # threshold is the average FDD of veqtl SNPs
p3 <- qtl_plot("eqtl_trans/new/MAF_map_noneqtl/Ctrl_FDD_summary.csv", "Ctrl trans-eQTL", 0.5212883999999999,"#4F8E4D")# threshold is the average FDD of eqtl SNPs
p4 <- qtl_plot("eqtl_trans/new/MAF_map_noneqtl/HS_FDD_summary.csv", "HS trans-eQTL", 0.5196066,"#4F8E4D")# threshold is the average FDD of eqtl SNPs

ggsave("MAF_veQTL_Ctrl_FDD_plot.svg", p1, width = 4, height = 4,dpi=300)
ggsave("MAF_veQTL_HS_FDD_plot.svg", p2, width = 4, height = 4,dpi=300)
ggsave("MAF_eQTL_Ctrl_FDD_plot.svg", p3, width = 4, height = 4,dpi=300)
ggsave("MAF_eQTL_HS_FDD_plot.svg", p4, width = 4, height = 4,dpi=300)

# cis
p5 <- qtl_plot("eqtl_cis/MAF_map/Ctrl_FDD_summary.csv", "Ctrl cis-eQTL", 0.49078271965989545, "#4F8E4D") # threshold is the FDD without between
p6 <- qtl_plot("eqtl_cis/MAF_map/HS_FDD_summary.csv", "HS cis-eQTL", 0.4929598408973634, "#4F8E4D") # threshold is the FDD without between
ggsave("eQTL_cis_Ctrl_FDD_plot.svg", p5, width = 4, height = 4,dpi=300)
ggsave("eQTL_cis_HS_FDD_plot.svg", p6, width = 4, height = 4,dpi=300)
"""
 
