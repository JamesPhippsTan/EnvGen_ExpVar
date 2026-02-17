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
