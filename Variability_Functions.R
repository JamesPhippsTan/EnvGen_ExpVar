# Variability functions
# A suite of functions for carrying out variability analyses

# Packages
library(matrixStats)
library(tidyverse)
library(scran) 
library(bootstrap) 
library(gridExtra)
library(dplyr)
library(ggplot2)
library(ggridges)


#### Rounding statistics on plots #####

format_statistic <- function(stat, threshold = 0.001) {
  if (is.na(stat)) {
    return("NA")
  } else if (abs(stat) == 0) { # some p < 2.2e-16 will store as 0
    return('< 2.2e-16')  
  } else if (abs(stat) < threshold) { # if below threshold, e format with 2 sig. figs
    return(formatC(stat, format = "e", digits = 2))
  } else {
    return(round(stat, 3)) # otherwise, 3 digits
  }
}

#### Plot overlapping proportion histograms for two or more different sets of genes ####
# Useful, for example, if one wants to plot the distribution of HS logCPMs of a subset of genes to that of all genes
# Takes the MeanVar dataframe, the gene set of interest and the reference set
# Outputs 2 histograms overlaid on one another
# Original goal: see the distribution of logCPMs of DVGs relative to the distribution of ALL genes

plot_twoset_gene_ridges <- function(df, gene_list1, gene_list2, gene_column, value_column) {
  # Subset and tag groups
  df1 <- df %>% 
    filter(!!sym(gene_column) %in% gene_list1) %>%
    mutate(group = "Set 1")
  
  df2 <- df %>% 
    filter(!!sym(gene_column) %in% gene_list2) %>%
    mutate(group = "Set 2")
  
  # Combine
  df_combined <- bind_rows(df1, df2)
  
  # Plot
  two_gene_set_prop_plot <- ggplot(df_combined, aes_string(x = value_column, fill = "group")) +
    geom_histogram(aes(y = ..count../sum(..count..)), 
                   position = "identity", alpha = 0.5, bins = 10) +
    labs(y = "Proportion of all genes", x = value_column, fill = "Group") +
    theme_classic()
  
  # Return
  (two_gene_set_prop_plot)
}

plot_multiset_gene_ridges <- function(
    df,
    gene_sets,        # named list: list(Set_A = c("gene1", "gene2"), Set_B = c("gene3", ...))
    gene_column,
    value_column,
    x_label = NULL,
    set_label = "Gene Set",
    palette = NULL,
    collapse_label=NULL
) {
  # Prepare combined data with group labels (replace _ with space in set names)
  data_list <- lapply(names(gene_sets), function(set_name) {
    display_name <- gsub("_", " ", set_name)
    if (!is.null(collapse_label) && collapse_label) {
      # Insert line breaks after specific words (customize as needed)
      display_name <- gsub("mean", "mean\n", display_name)
      display_name <- gsub("variability", "variability\n", display_name)
      display_name <- gsub("genes", "genes\n", display_name)
    } 
    genes <- gene_sets[[set_name]]
    df %>%
      filter(.data[[gene_column]] %in% genes) %>%
      mutate(group = paste0(display_name, " (n = ", n(), ")"))
  })
  df_combined <- bind_rows(data_list)
  
  # Set color palette if not provided
  if (is.null(palette)) {
    palette <- scales::hue_pal()(length(gene_sets))
    names(palette) <- unique(df_combined$group)
  }
  
  # Plot
  ggplot(df_combined, aes(
    x = .data[[value_column]],
    y = group,
    fill = group
  )) +
    geom_density_ridges(
      alpha = 0.8,
      scale = 1.5,
      rel_min_height = 0.0005,
      bandwidth = 0.6
    ) +
    scale_fill_manual(values = palette) +
    labs(
      x = if (!is.null(x_label)) x_label else value_column,
      y = NULL,
      fill = set_label
    ) +
    theme_classic() +
    theme(legend.position = "none")
}
#### Variability-mean relationship across the whole transcriptome ####

# Used to visualize the transcriptome-wide relationship between variability and mean
# In the form of a var~mean plot and rank(var)~rank(mean) plot
# For a subset of samples (e.g., condition-specific relationship)
# Variability is by default variance but can be changed to MAD
# Mean is by default mean but can be changed to median
# Log options available for visualisation

Variability_mean_plots <- function(exp_matrix,
                                   samples,
                                   var_func,
                                   mean_func,
                                   log_mean,
                                   log_var,
                                   mean_name,
                                   var_name){
  
  # Subset if needed (E.g., to view a condition-specific relationship)
  exp_matrix_samples <- as.matrix(exp_matrix[,samples])
  
  # Compute variability metric using the provided variability function
  Vars <- var_func(exp_matrix_samples)
  if (as.character(substitute(var_func)) == "rowMads"){
  Vars <- Vars/1.4826 # scale factor
  }

  # Compute mean/central tendency metric using the provided mean function
  Means <- mean_func(exp_matrix_samples)
  
  # Compute correlation tests
  # Pearson for actual values
  pe.test <- cor.test(y=Vars, x=Means, method = 'pearson')
  pe.coef <- pe.test$estimate
  pe.p.value <- pe.test$p.value
  # Spearman for ranks
  sp.test <- cor.test(y=Vars, x=Means, method = 'spearman')
  sp.coef <- sp.test$estimate
  sp.p.value <- sp.test$p.value

  # Make the two plots
  gg_df <- data.frame(Mean=Means,Var=Vars)
  
  # Var-mean rank plot
  pearson_stat_label <- paste0("Pearson corr. = ",format_statistic(pe.coef),
         ", p = ",format_statistic(pe.p.value))
  var_mean_rankplot <- ggplot(gg_df, aes(y = rank(Vars), x = rank(Means))) +
    geom_point(color = "#611BB8",alpha=0.5) +
    theme_classic() + 
    ggtitle(pearson_stat_label)+
    xlab(paste0(mean_name," (rank)")) + ylab(paste0(var_name," (rank)")) +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 7))
  # Var-mean plot (log if specified)
  if (log_var==T){
    Vars <- log(Vars)
  }
  if (log_mean==T){
    Means <- log(Means)
  }
  spearman_stat_label <- paste0("Spearman corr. = ",format_statistic(sp.coef),
         ", p = ",format_statistic(sp.p.value))
  var_mean_plot <- ggplot(gg_df, aes(y = Vars, x = Means)) +
    geom_point(color = "#611BB8",alpha=0.5)+
    theme_classic() + 
    ggtitle(spearman_stat_label)+
    xlab(mean_name) + ylab(var_name) +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 7)) 

  # Show both plots
  grid.arrange(var_mean_plot, var_mean_rankplot, ncol=2)
  
  # Return the individual plots as a list
  var_mean_plots <- list(var_mean_plot, var_mean_rankplot)
  return(var_mean_plots)
} 


# What if the means and variabilities are pre-computed in a different way?
Variability_mean_plots_precomp <- function(MeanVarTable,
                                   log_mean,
                                   log_var,
                                   mean_col,
                                   var_col,
                                   mean_name,
                                   var_name){
  
  # Obtain the variability and mean columns
  Vars <- MeanVarTable[,var_col]
  Means <-  MeanVarTable[,mean_col]
  
  # Compute correlation tests
  # Pearson for actual values
  pe.test <- cor.test(y=Vars, x=Means, method = 'pearson')
  pe.coef <- pe.test$estimate
  pe.p.value <- pe.test$p.value
  pearson_stat_label <- paste0("Pearson corr. = ",format_statistic(pe.coef),
                               ", p = ",format_statistic(pe.p.value))

  # Spearman for ranks
  sp.test <- cor.test(y=Vars, x=Means, method = 'spearman')
  sp.coef <- sp.test$estimate
  sp.p.value <- sp.test$p.value
  spearman_stat_label <- paste0("Spearman corr. = ",format_statistic(sp.coef),
                                ", p = ",format_statistic(sp.p.value))
  
  
  # Make the two plots
  gg_df <- data.frame(Mean=Means,Var=Vars)
  
  # Var-mean rank plot
  var_mean_rankplot <- ggplot(gg_df, aes(y = rank(Vars), x = rank(Means))) +
    geom_point(color = "#611BB8",alpha=0.5) +
    theme_classic() + 
    ggtitle(pearson_stat_label)+
    xlab(paste0(mean_name," (rank)")) + ylab(paste0(var_name," (rank)")) +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 7))
  # Var-mean plot (log if specified)
  if (log_var==T){
    Vars <- log(Vars)
  }
  if (log_mean==T){
    Means <- log(Means)
  }
  var_mean_plot <- ggplot(gg_df, aes(y = Vars, x = Means)) +
    geom_point(color = "#611BB8",alpha=0.5)+
    theme_classic() + 
    ggtitle(spearman_stat_label)+
    xlab(mean_name) + ylab(var_name) +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 7)) 
  
  # Show both plots
  grid.arrange(var_mean_plot, var_mean_rankplot, ncol=2)
  
  # Return the individual plots as a list
  var_mean_plots <- list(var_mean_plot, var_mean_rankplot)
  return(var_mean_plots)
} 


# Find where the subset of genes lies within the variance-mean space
subset_overlay_plots <- function(df,mean_col,var_col,subset_of_genes,subset_name){
  intersecting_genes <- intersect(subset_of_genes,df$gene_ID)
  rownames(df) <- df$gene_ID
  df <- df[,c(mean_col,var_col)]
  par(mfcol=c(1,2))
  title <- paste0(as.character(subset_name)," (",length(intersecting_genes)," Genes)")
  plot(y=df[,var_col],x=df[,mean_col],cex=2,ylab=var_col,xlab=mean_col,col='black',ylim=c(min(df[,var_col]),max(df[,var_col])),xlim=c(min(df[,mean_col]),max(df[,mean_col])),main=title)
  par(new = TRUE)
  plot(y=df[intersecting_genes,var_col],x=df[intersecting_genes,mean_col],cex=2,pch = 16,col='red',ylab=var_col,xlab=mean_col,ylim=c(min(df[,var_col]),max(df[,var_col])),xlim=c(min(df[,mean_col]),max(df[,mean_col])))
  df$var_rank <- rank(df[,var_col])
  df$mean_rank <- rank(df[,mean_col])
  plot(y=df[,"var_rank"],x=df[,"mean_rank"],cex=2,col='black',ylab=paste0(var_col,' Rank'),xlab=paste0(mean_col," Rank"),ylim=c(min(df[,"var_rank"]),max(df[,"var_rank"])),xlim=c(min(df[,"mean_rank"]),max(df[,"mean_rank"])))
  par(new = TRUE)
  plot(y=df[intersecting_genes,"var_rank"],x=df[intersecting_genes,"mean_rank"],cex=2,pch = 16,col='red',ylab=paste0(var_col,' Rank'),xlab=paste0(mean_col," Rank"),ylim=c(min(df[,"var_rank"]),max(df[,"var_rank"])),xlim=c(min(df[,"mean_rank"]),max(df[,"mean_rank"])))
}

subset_overlay_plots_object <- function(df, mean_col, var_col, subset_of_genes) {
  intersecting_genes <- intersect(subset_of_genes, df$gene_ID)
  rownames(df) <- df$gene_ID
  df <- df[, c(mean_col, var_col)]
  
  # Title for the plots
  title <- paste0(as.character(substitute(subset_of_genes)), " (", length(intersecting_genes), " Genes)")
  
  # Create an empty list to store plots
  plot_list <- list()
  
  # First plot: Raw values
  par(mfcol = c(1, 2))  # 1 row, 2 columns layout
  plot(y = df[, var_col], x = df[, mean_col], ylab = var_col, xlab = mean_col, col = 'black',
       ylim = range(df[, var_col]), xlim = range(df[, mean_col]), main = title)
  par(new = TRUE)
  plot(y = df[intersecting_genes, var_col], x = df[intersecting_genes, mean_col], col = 'red', pch = 16,
       ylab = "", xlab = "", axes = FALSE, ylim = range(df[, var_col]), xlim = range(df[, mean_col]))
  
  plot_list$raw_plot <- recordPlot()  # Save the first plot
  
  # Compute rank transformation
  df$var_rank <- rank(df[, var_col])
  df$mean_rank <- rank(df[, mean_col])
  
  # Second plot: Rank-transformed values
  par(mfcol = c(1, 2))
  plot(y = df[, "var_rank"], x = df[, "mean_rank"], col = 'black',
       ylab = paste0(var_col, ' Rank'), xlab = paste0(mean_col, " Rank"),
       ylim = range(df[, "var_rank"]), xlim = range(df[, "mean_rank"]))
  par(new = TRUE)
  plot(y = df[intersecting_genes, "var_rank"], x = df[intersecting_genes, "mean_rank"], col = 'red', pch = 16,
       ylab = "", xlab = "", axes = FALSE, ylim = range(df[, "var_rank"]), xlim = range(df[, "mean_rank"]))
  
  plot_list$ranked_plot <- recordPlot()  # Save the second plot
  
  return(plot_list)
}

#### Special variability metrics ####

# Adjusted variability from Liu et al. 2020 BMC Biology, inspired by Barosso et al., 2018 Genetics
# Based on a polynomial regression of variability against mean
polynomial_adjust_var <-function(global_means,global_vars,specific_vars) {
  expData <- data.frame(mean=global_means,var=global_vars,specific_vars=specific_vars)
  m1 <- lm(var~mean, expData)
  m2 <- update(m1, .~. + I(mean^2), expData)
  m3 <- update(m2, .~. + I(mean^3), expData)
  m4 <- update(m3, .~. + I(mean^4), expData)
  m5 <- update(m4, .~. + I(mean^5), expData)
  m6 <- update(m5, .~. + I(mean^6), expData)
  m7 <- update(m6, .~. + I(mean^7), expData)
  m8 <- update(m7, .~. + I(mean^8), expData)
  m9 <- update(m8, .~. + I(mean^9), expData)
  m10 <- update(m9, .~. + I(mean^10), expData)
  totalM<-list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
  for (i in c(1:9)) {
    modelComp <- anova(totalM[[i]],totalM[[i+1]])
    if (is.na(modelComp$`Pr(>F)`[2]) | modelComp$`Pr(>F)`[2]>0.05) {
      m<-totalM[[i]]
      print(i)
      break
    }
  }
  adjusted_var <- expData$specific_vars/predict(m)
  return(list(Vars=adjusted_var,polynomial=i))
}

#### Jackknifing, a.k.a. leave-one-out estimation ####

# Make list of tables containing jackknife estimates and how far they are from the total population estimate
jackknife_tables <-function(exp_data,
                               samples,
                               func) {
  
  # Get the full population estimates
  exp_data <- as.matrix(exp_data[,samples])
  population_estimates <- func(exp_data)
  
  # Get the jackknifed estimates
  jackknife_estimates <- exp_data # retain the column and rownames
  for (sample in 1:ncol(exp_data)){
    jackknife_sample <- exp_data[,-sample] # remove one sample
    jackknife_estimates[,sample] <- func(as.matrix(jackknife_sample)) # calculate the new estimates
  }
  
  # If the function is MAD, divide by a constant
  if (as.character(substitute(var_func)) == "rowMads"){
    population_estimates <- population_estimate/1.4826
    jackknife_estimates <- jackknife_estimates/1.4826 
  }
  
  # Get the range of differences between the jackknifed estimates and full population estimate
  percent_differences  <- 100*(jackknife_estimates-population_estimates)/population_estimates
  percent_differences_range <- range(percent_differences)
  rank_differences  <- t(colRanks(jackknife_estimates))-rank(population_estimates)
  rank_differences_range <- range(rank_differences)
  
  # Save as a list object
  jackknife_table_list <- list(population_estimates=population_estimates,
                               jackknife_estimates=jackknife_estimates,
                               percent_differences=percent_differences,
                               percent_differences_range=percent_differences_range,
                               rank_differences=rank_differences,
                               rank_differences_range=rank_differences_range)
  return(jackknife_table_list)
}

