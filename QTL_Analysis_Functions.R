library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

# Function to split SNP name into chromosome and position
extract_chromosome_and_pos <- function(data, input_col, 
                                       chrom_col = "CHROM",
                                       pos_col = "POS",
                                       keep_original = FALSE) {
  # Define chromosome patterns (modify as needed)
  dros_chromosomes <- c("2L", "2R", "3L", "3R", "4", "X")
  pattern_regex <- str_c(dros_chromosomes, collapse = "|")
  
  # Process the data
  result <- data %>%
    mutate(
      # Extract chromosome (first match)
      !!chrom_col := str_extract(!!sym(input_col), pattern_regex),
      
      # Extract pos after chromosome
      !!pos_col := ifelse(
        is.na(!!sym(chrom_col)),
        NA_character_,
        str_remove(!!sym(input_col), str_c(".*?", !!sym(chrom_col)))
      )
    )
  # Remove original column if requested
  if (!keep_original) {
    result <- result %>% select(-all_of(input_col))
  }
  return(result)
}



# Functions
# Remove cis-veQTL from trans-veQTL mapping results
remove_cis_veQTL <- function(trans_mapping_df,SNP_positions,gene_positions,cis_window){
  # Add the actual SNP positions to the trans-mapping result
  trans_mapping_df <- merge(trans_mapping_df,SNPs_position_map,all.x=T,by="POS")
  # Add the positions of the gene bodies to the dataframe
  trans_mapping_df <- merge(trans_mapping_df,Gene_body_locations,all.x=T,by="GENE")
  # Call cis or trans
  trans_mapping_df$cis_or_trans <- ifelse(
    trans_mapping_df$GENE_CHR==trans_mapping_df$Actual_CHROM &
      trans_mapping_df$Actual_POS>(trans_mapping_df$GENE_START-cis_window) & 
      trans_mapping_df$Actual_POS<(trans_mapping_df$GENE_END+cis_window),
    'cis',
    'trans'
  )
  # Retain and return only trans-mapping tests
  trans_mapping_df_cis_removed <- subset(trans_mapping_df,cis_or_trans=='trans')
  return(trans_mapping_df_cis_removed)
}

# Map dummy SNP ID to actual SNP IDs based on a provided mapping between them
# Sometimes the SNP IDs vary depending on the dataset
get_real_SNP_ID <- function(dummy_IDs,mapping_df,dummy_ID_column,real_ID_column){
  mapping_df <- column_to_rownames(mapping_df,var=dummy_ID_column)
  real_IDs <- mapping_df[dummy_IDs,real_ID_column]
  return(real_IDs)
}

# Number of veQTLs per gene
veQTL_number_per_gene <- function(gene_df){
  split_gene_df <- split(gene_df,f=gene_df[,'GENE'])
  number_of_veQTL_by_gene <- as.data.frame(do.call(rbind,lapply(split_gene_df, nrow)))
  number_of_veQTL_by_gene$GENE <- rownames(number_of_veQTL_by_gene)
  return(number_of_veQTL_by_gene)
}

# Create tables with summaries of the number of SNPs or genes per hotspot genes and hotspot SNPs
summary_with_count <- function(col) {
  sum_stats <- summary(col)
  num_non_na <- sum(!is.na(col))
  c(sum_stats, Total_number = num_non_na)
}

# Handle the lower limit of p-values
format_P_handle0 <- function(p) {
  if (is.na(p)) {
    return("NA")
  } else if (p == 0) {
    return("2.2e-16")
  } else {
    return(formatC(p, format = "e", digits = 2))
  }
}

# Add SNP stats and effects to an existing QTL dataframe
add_SNP_stats <- function(QTL_df,SNP_stats_df,merge_column){
  SNP_stats_df$variant_id <- SNP_stats_df$SNP
  # Get sample size
  sample_size <- SNP_stats_df[,"C(HOMA1)"] + SNP_stats_df[,"C(HOMA2)"] + SNP_stats_df[,"C(HET)"]
  # Calculate MAF
  allele_freq <- (2 * SNP_stats_df[,"C(HOMA1)"] + SNP_stats_df[,"C(HET)"]) / (2 * sample_size)
  SNP_stats_df$MAF <- pmin(allele_freq, 1 - allele_freq)
  # Calculate genotype frequencies
  SNP_stats_df$ALT_freq <- SNP_stats_df[,"C(HOMA1)"]/sample_size
  SNP_stats_df$REF_freq <- SNP_stats_df[,"C(HOMA2)"]/sample_size
  SNP_stats_df$HET_freq <- SNP_stats_df[,"C(HET)"]/sample_size
  # Calculate sample size and missingness
  SNP_stats_df$Sample_size <- sample_size
  SNP_stats_df$Missingness <- SNP_stats_df[,"C(MISSING)"]/(SNP_stats_df[,"C(MISSING)"]+sample_size)
  # Minor_homozygote_count
  SNP_stats_df$Minor_hom_count <- pmin(SNP_stats_df[,"C(HOMA1)"], SNP_stats_df[,"C(HOMA2)"])
  SNP_stats_df$Minor_hom_het_count <- pmin(SNP_stats_df[,"C(HOMA1)"], SNP_stats_df[,"C(HOMA2)"])+SNP_stats_df[,"C(HET)"]
  # Merge QTL df and stats
  QTL_df <- merge(x=QTL_df,y=SNP_stats_df,all.x=T,by=merge_column)
  # Merge
  return(QTL_df)
}

# Function to create overlap plots (horizontal bars) as an alternative to Venn diagrams
# Requires a dataframe with the overlaps and non-overlap numbers pre-calculated 
# Can be used for any feature, e.g., eGenes, vGenes, veQTL, eQTL, etc.
create_overlap_plot <- function(data,partition_colors,feature) {
  ggplot(data, aes(x = Value, y = Partition, fill = Partition)) +
    geom_col(width = 0.6) +
    geom_text(
      aes(label = paste0(round(Percent * 100, digits = 1), "%")),
      hjust = -0.1,
      size = 4
    ) +
    scale_fill_manual(values = partition_colors) +
    labs(
      x = NULL,
      y = NULL,
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      axis.title.x.top = element_text(),
      axis.text.x.top = element_text(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    scale_x_continuous(
      name = feature,
      position = "top",
      expand = expansion(mult = c(0, 0.4)))
}

# Add MAF
add_MAF <- function(SNP_stats_df){
  sample_size <- SNP_stats_df[,"C(HOMA1)"] + SNP_stats_df[,"C(HOMA2)"] + SNP_stats_df[,"C(HET)"]
  # Calculate MAF
  # Frequency of Allele 1 (A1)
  allele_freq <- (2 * SNP_stats_df[,"C(HOMA1)"] + SNP_stats_df[,"C(HET)"]) / (2 * sample_size)
  # Minor allele frequency (regardless if A1 is minor or not)
  SNP_stats_df$MAF <- pmin(allele_freq, 1 - allele_freq)
  # Determine if A1 is or is not the minor allele 
  SNP_stats_df$Minor_allele <- ifelse(allele_freq<0.5,SNP_stats_df$A1,SNP_stats_df$A2)
  # If truly 50:50, then cannot call
  SNP_stats_df$Minor_allele <- ifelse(allele_freq==0.5,'NA',SNP_stats_df$Minor_allele)
  return(SNP_stats_df)
}

# Plot eQTL and veQTL MAF histograms next to each other
compare_MAFs_eQTL_veQTL <- function(eQTL_datasets, veQTL_datasets, MAF_df) {
  # MAF column selector
  get_maf_column <- function(name) {
    if (grepl("^Ctrl", name)) return("Ctrl_MAF")
    if (grepl("^HS", name)) return("HS_MAF")
    return("Avg_MAF")  # For "All" sets
  }
  
  results <- data.frame(
    Comparison = character(),
    KS_p = numeric(),
    T_p = numeric(),
    medianDiff = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(eQTL_datasets)) {
    if (!name %in% names(veQTL_datasets)) next
    
    message("Processing: ", name)
    maf_col <- get_maf_column(name)
    
    eqtl_snps <- unique(eQTL_datasets[[name]])
    veqtl_raw <- veQTL_datasets[[name]]
    veqtl_snps <- if (is.data.frame(veqtl_raw)) unique(veqtl_raw$SNP) else unique(veqtl_raw)
    
    # Filter to SNPs present in MAF_df
    eqtl_snps <- eqtl_snps[eqtl_snps %in% rownames(MAF_df)]
    veqtl_snps <- veqtl_snps[veqtl_snps %in% rownames(MAF_df)]
    
    eqtl_maf <- MAF_df[eqtl_snps, maf_col, drop = TRUE]
    veqtl_maf <- MAF_df[veqtl_snps, maf_col, drop = TRUE]
    
    # Remove NA values
    eqtl_maf <- eqtl_maf[!is.na(eqtl_maf)]
    veqtl_maf <- veqtl_maf[!is.na(veqtl_maf)]
    
    if (length(eqtl_maf) < 3 || length(veqtl_maf) < 3) next
    
    # Statistical tests
    ks_res <- ks.test(eqtl_maf, veqtl_maf)
    w_res <- wilcox.test(eqtl_maf, veqtl_maf)
    median_eqtl <- median(eqtl_maf)
    median_veqtl <- median(veqtl_maf)
    
    # Format p-values
    format_P_handle0 <- function(p) {
      if (is.na(p)) return("NA")
      else if (p == 0) return("2.2e-16")
      else return(formatC(p, format = "e", digits = 2))
    }
    ks_p <- format_P_handle0(ks_res$p.value)
    wilcox_p <- format_P_handle0(w_res$p.value)
    
    results <- rbind(results, data.frame(
      Comparison = name,
      KS_p = ks_p,
      Wilcox_p = wilcox_p,
      eQTL_median_MAF = median_eqtl,
      veQTL_median_MAF = median_veqtl,
      eQTL_n = length(eqtl_maf),
      veQTL_n = length(veqtl_maf)
    ))
    
    # Prepare plot data
    plot_data <- rbind(
      data.frame(Group = "eQTL", MAF = eqtl_maf),
      data.frame(Group = "veQTL", MAF = veqtl_maf)
    )
    
    # Bin and normalize within each group
    bin_breaks <- seq(0, 0.5, length.out = 11)
    plot_data_binned <- plot_data %>%
      mutate(Bin = cut(MAF, breaks = bin_breaks, include.lowest = TRUE, right = FALSE)) %>%
      group_by(Group, Bin) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Group) %>%
      mutate(Percent = Count / sum(Count) * 100)
    
    # Plot with all requested changes
    p_legend <- ggplot(plot_data_binned, aes(x = Bin, y = Percent, fill = Group)) +
      geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
      theme_classic() +
      labs(
        title = paste0("Kolmogorov-Smirnoff p = ", ks_p, " | Wilcox p = ", wilcox_p),
        x = "MAF",
        y = "Percent of all QTL"
      ) +
      scale_fill_manual(
        values = c("eQTL" = "#4F8E4D", "veQTL" = "#611BB8"),  # Updated eQTL color
        labels = c(
          paste0("eQTL (n = ", length(eqtl_maf), ", median = ", round(median_eqtl, 3), ")"),
          paste0("veQTL (n = ", length(veqtl_maf), ", median = ", round(median_veqtl, 3), ")")
        ),
        name = NULL
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = c(1, 1),  # Legend inside plot (top-right)
        legend.justification = c(1, 1),  # Anchor to top-right corner
        legend.background = element_rect(fill = "white", color = "black"),  # Legible background
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
      )
    p <- p_legend + theme(legend.position = "none")
    setwd("C:\\Users\\jtanshengyi\\Desktop\\Projects\\veQTL Netherlands Normal vs High Sugar Adult\\Data\\GraVe_Mapping\\MAF")
    ggsave(filename = paste0(name, "_MAF_comparison.svg"), plot = p, width = 4.7, height = 4,dpi=300)
    ggsave(filename = paste0(name, "_MAF_comparison_legend.svg"), plot = p_legend, width = 4.7, height = 4,dpi=300)
  }
  return(results)
}

########################### Regulatory elements #############################

# Function to process all the scEnhancer text files within a directory
# into one dataframe with the relevant columns for downstream analysis
process_scEnhancer_tables <- function(directory) {
  # Get all .txt files in the directory
  file_list <- list.files(path = directory, pattern = "\\.txt$", full.names = TRUE)
  
  # Read each file and add filename as a column
  enhancer_list <- lapply(file_list, function(file) {
    df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    df$source_file <- basename(file)
    return(df)
  })
  
  # Combine all data frames into one
  scEnhancers <- do.call(rbind, enhancer_list)
  
  # Rename columns
  colnames(scEnhancers)[1:2] <- c("enhancer_info", "gene_info")
  
  # Split enhancer_info column (more stable version)
  scEnhancers <- scEnhancers %>%
    tidyr::extract(
      col = enhancer_info,
      into = c("CHROM", "START", "END", "enhancer_score"),
      regex = "^([^:]+):([^-]+)-([^|]+)\\|([^|]+)",
      remove = TRUE
    )
  
  # Split gene_info column with more robust handling
  scEnhancers <- scEnhancers %>%
    tidyr::extract(
      col = gene_info,
      into = c("gene_CHROM", "gene_START", "gene_END", "gene_identifiers"),
      regex = "^([^:]+):([^-]+)-([^|]+)\\|([^|]*)",
      remove = TRUE
    ) %>%
    dplyr::mutate(
      # Handle gene identifiers more carefully
      gene_identifier = sapply(strsplit(gene_identifiers, ":"), function(x) ifelse(length(x) > 0, x[1], NA)),
      gene_identifier2 = sapply(strsplit(gene_identifiers, ":"), function(x) ifelse(length(x) > 1, x[2], NA)),
      # Clean up gene_identifier2 (remove trailing | if present)
      gene_identifier2 = sub("\\|$", "", gene_identifier2)
    ) %>%
    dplyr::select(-gene_identifiers)  # Remove the temporary column
  
  # Convert numeric columns to appropriate types
  scEnhancers <- scEnhancers %>%
    dplyr::mutate(
      START = as.integer(START),
      END = as.integer(END),
      enhancer_score = as.numeric(enhancer_score),
      gene_START = as.integer(gene_START),
      gene_END = as.integer(gene_END)
    )
  
  # Reorder columns
  final_columns <- c("CHROM", "START", "END", "enhancer_score",
                     "gene_CHROM", "gene_START", "gene_END",
                     "gene_identifier", "gene_identifier2",
                     "source_file")
  
  scEnhancers <- scEnhancers %>%
    dplyr::select(dplyr::any_of(final_columns), dplyr::everything())
  
  scEnhancers$CHROM <- str_remove(scEnhancers$CHROM, "^chr")
  scEnhancers$CHROM_wHet <- scEnhancers$CHROM
  scEnhancers$CHROM <- str_remove(scEnhancers$CHROM, "Het")
  scEnhancers$name <-paste0(gsub('interaction.txt','',scEnhancers$source_file),scEnhancers$CHROM_wHet,scEnhancers$START)
  return(scEnhancers)
}

# Get all SNPs per category (excluding NA and 0)
create_category_snps <- function(snp_df) {
  # Exclude first column (e.g., "ID") and SNP column
  category_names <- names(snp_df)[-c(1, 2)]  # Adjust if SNP is not column 2
  snp_lists <- list()
  
  for (category in category_names) {
    # Get SNPs with non-NA and non-zero values
    snps <- snp_df %>%
      filter(!is.na(.data[[category]]), .data[[category]] != 0) %>%
      pull(SNP) %>%
      unique()
    
    snp_lists[[category]] <- snps
  }
  return(snp_lists)
}

# Get top N SNPs per category. 
# Can either be unique (have eGenes/vGenes only in the condition) or 
# shared (>1 regulated genes in two or more conditions)
create_top_snps <- function(snp_df, 
                            n = 10, 
                            mode = c("shared", "unique")) {
  # Validate inputs
  mode <- match.arg(mode)
  category_names <- names(snp_df)[-1]  # Exclude SNP column
  
  # Calculate sum of other categories (for "shared" mode only)
  snp_df <- snp_df %>%
    mutate(
      other_categories_sum = if (mode == "shared") {
        rowSums(select(., -1), na.rm = TRUE)
      } else {
        NA_real_  # Not used in "unique" mode
      }
    )
  
  top_snp_lists <- list()
  
  for (category in category_names) {
    if (mode == "unique") {
      # Get SNPs where:
      # 1. Current category is non-NA/non-zero
      # 2. All other categories are NA
      candidate_snps <- snp_df %>%
        filter(
          !is.na(.data[[category]]),
          .data[[category]] != 0,
          if_all(-1 & !all_of(category), ~is.na(.))  # All other columns must be NA
        )
    } else {
      candidate_snps <- snp_df %>%
        filter(!is.na(.data[[category]]), .data[[category]] != 0)
    }
    
    if (nrow(candidate_snps) == 0) {
      top_snp_lists[[category]] <- character(0)
      next
    }
    
    if (mode == "shared") {
      ranked_snps <- candidate_snps %>%
        arrange(desc(.data[[category]]), desc(other_categories_sum))
      
      if (nrow(ranked_snps) >= n) {
        nth_value <- ranked_snps[[category]][n]
        nth_sum <- ranked_snps$other_categories_sum[n]
        top_snps <- ranked_snps %>%
          filter(.data[[category]] > nth_value | 
                   (.data[[category]] == nth_value & other_categories_sum >= nth_sum))
      } else {
        top_snps <- ranked_snps
      }
    } else {
      top_snps <- candidate_snps %>%
        arrange(desc(.data[[category]])) %>%
        slice_head(n = n)
    }
    
    top_snp_lists[[category]] <- top_snps$SNP
  }
  
  return(top_snp_lists)
}

# Function to search a region around the TSS of a promoter for the SNP of interest
promoter_search <- function(snp_list, promoter_df, window_size = 50) {
  # Initialize results data frame
  results <- data.frame(
    snp_name = character(),
    chromosome = character(),
    snp_position = numeric(),
    tss_position = numeric(),
    difference = numeric(),
    gene_name = character(),
    near_tss = character(),
    stringsAsFactors = FALSE
  )
  
  # Define chromosome patterns (2L, 2R, 3L, 3R, 4, X)
  chrom_patterns <- c("2L", "2R", "3L", "3R", "4", "X")
  
  # Process each SNP
  for (snp in snp_list) {
    # Extract chromosome and position from SNP name
    chrom <- NULL
    snp_pos <- NA
    
    # Try each chromosome pattern to find a match
    for (pattern in chrom_patterns) {
      if (grepl(paste0("^", pattern), snp)) {
        chrom <- pattern
        snp_pos <- as.numeric(sub(paste0("^", pattern), "", snp))
        break
      }
    }
    
    if (is.null(chrom)) {
      warning(paste("SNP", snp, "doesn't match expected chromosome patterns"))
      next
    }
    
    if (is.na(snp_pos)) {
      warning(paste("SNP", snp, "doesn't contain a valid position"))
      next
    }
    
    # Get TSS positions for this chromosome
    chrom_tss <- promoter_df[promoter_df$CHROM == chrom, ]
    
    if (nrow(chrom_tss) == 0) {
      # No TSS for this chromosome
      results <- rbind(results, data.frame(
        snp_name = snp,
        chromosome = chrom,
        snp_position = snp_pos,
        tss_position = NA,
        difference = NA,
        gene_name = NA,
        near_tss = "no",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Initialize variables to track the closest TSS
    closest_diff <- Inf
    closest_tss <- NA
    closest_gene <- NA
    near_tss_flag <- FALSE
    
    # Check each TSS for this chromosome
    for (i in 1:nrow(chrom_tss)) {
      current_tss <- chrom_tss$START[i]
      current_gene <- chrom_tss$gene_name[i]
      diff <- snp_pos - current_tss
      abs_diff <- abs(diff)
      
      # Update closest TSS if this one is closer
      if (abs_diff < abs(closest_diff)) {
        closest_diff <- diff
        closest_tss <- current_tss
        closest_gene <- current_gene
      }
      
      # Check if within window
      if (abs_diff <= window_size) {
        near_tss_flag <- TRUE
      }
    }
    
    # Add to results
    results <- rbind(results, data.frame(
      snp_name = snp,
      chromosome = chrom,
      snp_position = snp_pos,
      tss_position = closest_tss,
      difference = closest_diff,
      gene_name = closest_gene,
      near_tss = ifelse(near_tss_flag, "yes", "no"),
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

genomic_element_search <- function(snp_list, genomic_element_df) {
  # Initialize results data frame
  results <- data.frame(
    snp_name = character(),
    chromosome = character(),
    snp_position = numeric(),
    element_name = character(),
    gene_identifier = character(),
    overlaps = character(),
    stringsAsFactors = FALSE
  )
  
  # Define chromosome patterns (2L, 2R, 3L, 3R, 4, X)
  chrom_patterns <- c("2L", "2R", "3L", "3R", "4", "X")
  
  # Process each SNP
  for (snp in snp_list) {
    # Extract chromosome and position from SNP name
    chrom <- NULL
    snp_pos <- NA
    
    # Try each chromosome pattern to find a match
    for (pattern in chrom_patterns) {
      if (grepl(paste0("^", pattern), snp)) {
        chrom <- pattern
        snp_pos <- as.numeric(sub(paste0("^", pattern), "", snp))
        break
      }
    }
    
    if (is.null(chrom)) {
      warning(paste("SNP", snp, "doesn't match expected chromosome patterns"))
      next
    }
    
    if (is.na(snp_pos)) {
      warning(paste("SNP", snp, "doesn't contain a valid position"))
      next
    }
    
    # Get genomic elements for this chromosome
    chrom_elements <- genomic_element_df[genomic_element_df$CHROM == chrom, ]
    
    if (nrow(chrom_elements) == 0) {
      # No elements for this chromosome
      results <- rbind(results, data.frame(
        snp_name = snp,
        chromosome = chrom,
        snp_position = snp_pos,
        element_name = NA,
        gene_identifier = NA,
        overlaps = "no",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Check if SNP falls within any element's range
    overlapping_elements <- chrom_elements[
      chrom_elements$START <= snp_pos & chrom_elements$END >= snp_pos,
    ]
    
    if (nrow(overlapping_elements) > 0) {
      # SNP overlaps with at least one element
      # We'll take the first match (or you could modify to handle multiple matches)
      results <- rbind(results, data.frame(
        snp_name = snp,
        chromosome = chrom,
        snp_position = snp_pos,
        element_name = overlapping_elements[,'element_name'][1],
        gene_identifier = overlapping_elements[,'gene_identifier'][1],
        overlaps = "yes",
        stringsAsFactors = FALSE
      ))
    } else {
      # No overlap found
      results <- rbind(results, data.frame(
        snp_name = snp,
        chromosome = chrom,
        snp_position = snp_pos,
        element_name = NA,
        gene_identifier = NA,
        overlaps = "no",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(results)
}


# Function apply all four types of search to a list of SNP sets
process_snp_sets <- function(snp_sets, analysis_functions) {
  all_results <- list()
  feature_fractions <- list()
  
  for (snp_set_name in names(snp_sets)) {
    snp_list <- snp_sets[[snp_set_name]]
    set_results <- list()
    set_fractions <- list()
    
    for (func_name in names(analysis_functions)) {
      # Apply the current function
      result <- analysis_functions[[func_name]](snp_list)
      
      # Store the complete results
      result_name <- paste(snp_set_name, func_name, sep = "_")
      all_results[[result_name]] <- result
      
      # Calculate the fraction with feature overlap
      if (func_name == "promoter_proximity") {
        frac <- mean(result$near_tss == "yes", na.rm = TRUE)
      } else {
        frac <- mean(result$overlaps == "yes", na.rm = TRUE)
      }
      set_fractions[[func_name]] <- frac
    }
    
    # Store fractions for this SNP set
    feature_fractions[[snp_set_name]] <- data.frame(
      SNP_set = snp_set_name,
      as.data.frame(set_fractions)
    )
  }
  
  # Combine all fractions into one data frame
  fraction_df <- bind_rows(feature_fractions)
  
  return(list(
    detailed_results = all_results,
    feature_fractions = fraction_df
  ))
}


###################### SNP features enrichments using chi-square tests for multiple SNP sets ##############

# Functions to calculate how many yes and nos, as well as the ratio of yes to nos
get_counts <- function(df, type_col) {
  # Ensure the table includes both "yes" and "no" categories
  counts <- table(factor(df[[type_col]], levels = c("yes", "no")))
  return(counts)
}
calculate_yes_no_ratio <- function(counts) {
  ratio <- counts["yes"] / counts["no"]
  return(ratio)
}
# Create contingency tables comparing the number of members fo each category against all other background SNPs
create_contingency_table <- function(condition_SNPs, all_df, type_col) {
  # Get counts for the condition-specific data frame
  condition_counts <- get_counts(all_df[condition_SNPs,], type_col)
  
  # Get counts for the remaining SNPs (all SNPs minus condition-specific SNPs)
  all_counts <- get_counts(all_df, type_col)
  remaining_counts <- all_counts - condition_counts
  
  # Calculate ratios
  condition_ratio <- calculate_yes_no_ratio(condition_counts)
  remaining_ratio <- calculate_yes_no_ratio(remaining_counts)
  
  # Create the contingency table
  contingency_table <- rbind(condition_counts, remaining_counts)
  
  # Return the contingency table and ratios
  return(list(
    contingency_table = contingency_table,
    condition_counts = condition_counts,
    remaining_counts = remaining_counts,
    condition_ratio = condition_ratio,
    remaining_ratio = remaining_ratio
  ))}
# Then run chi-square tests
perform_chi_square_test <- function(contingency_table) {
  test_result <- chisq.test(contingency_table)
  return(test_result)
}


########################### Allele age #########################

# Fraction of SNPs where the derived allele increased transcript level or variability
Fraction_Derived_Increased <- function(QTL_df, SNP_ages_df, slope_allele_column, slope_column) {
  # Create unified SNP ID in SNP ages dataframe
  SNP_ages_df$SNP <- paste0(SNP_ages_df$arm, SNP_ages_df$snp_pos)
  
  # Merge with QTL data by SNP
  QTL_df <- merge(QTL_df, SNP_ages_df, by = "SNP", all = TRUE)
  
  # Define relevant columns
  alt_col <- slope_allele_column
  slope <- slope_column
  
  # Keep necessary columns only
  QTL_df <- QTL_df[, c("SNP", "ANC", alt_col, slope)]
  
  # Remove rows without ancestral allele assignment
  QTL_df <- QTL_df[!is.na(QTL_df[["ANC"]]), ]
  
  # Direction of effect with respect to derived allele
  QTL_df$DER_Direction <- ifelse(QTL_df[["ANC"]] != QTL_df[[alt_col]], 
                                 QTL_df[[slope]], 
                                 -QTL_df[[slope]])
  
  # Determine if derived allele increases trait
  QTL_df$Derived_Increased <- ifelse(QTL_df$DER_Direction > 0, 
                                     "Derived_Increased", 
                                     "Derived_Decreased")
  
  # Tabulate and compute percentage
  freq_table <- table(QTL_df$Derived_Increased)
  percentage_table <- (freq_table / sum(freq_table)) * 100
  
  # Ensure 0% is returned if no "Derived_Increased" is present
  return(if ("Derived_Increased" %in% names(percentage_table)) {
    percentage_table[["Derived_Increased"]]
  } else {
    message("Not a single slope in the dataframe is positive. Setting output to 0")
    0
  })
}

# Get SNP-specific FDI
get_SNP_specific_FDI <- function(QTL_df, SNP_ages_df, slope_allele_column, slope_column) {
  # Convert to data.table
  QTL <- as.data.table(QTL_df)
  SNP_ages <- as.data.table(SNP_ages_df)
  
  # Create SNP ID in SNP_ages_df
  SNP_ages[, SNP := paste0(arm, snp_pos)]
  
  # Merge to get ancestral state
  merged <- merge(QTL, SNP_ages[, .(SNP, ANC)], by = "SNP", all.x = TRUE)
  
  # Remove rows with missing ancestral state
  merged <- merged[!is.na(ANC)]
  
  # Calculate direction with respect to derived allele
  merged[, derived_direction := ifelse(get(slope_allele_column) != ANC, 
                                       get(slope_column), 
                                       -get(slope_column))]
  
  # Indicator for derived allele increasing the trait
  merged[, derived_increased := as.integer(derived_direction > 0)]
  
  # Summarise: number of genes and fraction derived increased (FDI)
  result <- merged[, .(
    Number_of_vGenes = .N,
    FDI = mean(derived_increased)
  ), by = SNP]
  
  return(result[])
}

# Get SNP-specific FMI
get_SNP_specific_FMI_HS <- function(QTL_df, Minor_allele_df, slope_allele_column, slope_column) {
  # Convert to data.table
  QTL <- as.data.table(QTL_df)
  SNP_minor_allele <- as.data.table(Minor_allele_df)

  # Merge to get minor allele 
  merged <- merge(QTL, SNP_minor_allele[, c("SNP","HS_minor_allele")], by = "SNP", all.x = TRUE)
  
  # Calculate direction with respect to minor allele
  merged[, minor_direction := ifelse(get(slope_allele_column) == HS_minor_allele, 
                                       get(slope_column), 
                                       -get(slope_column))]
  
  # Indicator for minor allele increasing the trait
  merged[, minor_increased := as.integer(minor_direction > 0)]
  
  # Summarise: number of genes and fraction minor increased (FMI)
  result <- merged[, .(
    Number_of_vGenes = .N,
    FMI = mean(minor_increased)
  ), by = SNP]
  
  return(result[])
}
