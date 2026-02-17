library(gprofiler2)
library(ggplot2)
library(rrvgo)
library(ggvenn)
library(gridExtra)
library(grid)
library(patchwork)

# Set the base URL - sometimes I cannot access the server if this is not set
set_base_url('https://biit.cs.ut.ee/gprofiler/')

##################################
##### gProfiler functions  #######
##################################


# GO enrichment analysis with custom background using gProfiler
perform_GO_enrichment <- function(gene_set, background_genes) {
  
  # Perform GO enrichment analysis with custom background gene set
  GO_result <- gost(query = gene_set, custom_bg = background_genes,
                    organism = "dmelanogaster",
                    evcodes = TRUE) # includes the intersecting genes
  GO_table <- data.frame(GO_result$result)
  # Return the result
  return(GO_table)
}

# Extract and format top GO terms (lowest FDR-corrected p-value) of each type
extract_top_GO_terms <- function(enrichment_result, top_n) {
  top_enriched_terms <- data.frame()
  for (GO_type in unique(enrichment_result$source)){
    result_df_GO_type <- subset(enrichment_result,source==GO_type)
    if (nrow(result_df_GO_type) >= top_n){
      result_df_GO_type <- result_df_GO_type[1:top_n,]
    }
    top_enriched_terms <- rbind(top_enriched_terms,result_df_GO_type)
  }
  return(top_enriched_terms)
}

# Generate gProfiler results for a list of gene sets
gProfiler_genesets <- function(gene_sets,background_genes) {
    
  all_results <- data.frame()
  top_20_results <- data.frame()
  
  # Loop through each gene set in the list
  for (dataset_name in names(gene_sets)) {
    if (is.null(rownames(gene_sets[[dataset_name]]))==T){ # If already a list of genes
      gene_set_genes <- gene_sets[[dataset_name]]
    }
    else{
      gene_set_genes <- rownames(gene_sets[[dataset_name]]) # If embedded in rownames of a file
    }
    number_of_genes <- length(gene_set_genes)
    print(paste0('Dataset ',dataset_name,' has ',as.character(number_of_genes),' Genes'))
    
    if (number_of_genes<=1){
      print('Fewer than 2 genes')
      next
    }
    
    result <- tryCatch(perform_GO_enrichment(gene_set_genes, background_genes = background_genes) , error= function(e) {return(0)}  )
    
    if (ncol(result)==0) {
      print("No enrichments")
      next
    }
    
    # Add dataset name to the dataframe
    result$dataset <- dataset_name
    
    # Append the results to the existing dataframe
    all_results <- rbind(all_results, result)
    
    # Top 20 GO terms of interest
    top_GO_terms <- extract_top_GO_terms(enrichment_result=result,top_n=20)
    top_20_results <- rbind(top_20_results, top_GO_terms)
    }
  return(list(all_results,top_20_results))
}

# Plot significant GO terms for a specific GO type
subset_GO = function (data, GO_type, dataset_name) {
  dataset <- data[[dataset_name]]
  dataset <- subset(dataset,source==GO_type)
  return(dataset)
}

# Plot significant top 20 significant GO terms 
plot_top_GO = function (data, GO_type, dataset_name) {
  dataset <- data[[dataset_name]]
  dataset <- subset(dataset,source==GO_type)
  ggplot(dataset) +
    geom_point(aes(x = (intersection_size/term_size)/(term_size/number_of_expressed_genes), y = term_name, 
                   color = p_value))+ 
    labs(y = paste0("Top 20 ",GO_type," terms"), x = "Gene Ratio",color= "p.adjust")+
    ggtitle(dataset_name)
}

# Refined version of plot_top_GO
plot_top_GO_refined = function(data, GO_type, dataset_name, 
                       x_lab = "Gene Ratio", 
                       max_label_length = 20, 
                       max_lines = 2) {
  
  # Check if dataset exists
  if (!dataset_name %in% names(data)) {
    stop("Dataset '", dataset_name, "' not found in the provided data")
  }
  
  # Extract and filter dataset
  dataset <- data[[dataset_name]]
  dataset <- subset(dataset, source == GO_type)

  # Wrap long term names into multiple lines
  wrap_labels <- function(labels, length_limit = max_label_length, line_limit = max_lines) {
    sapply(labels, function(label) {
      if (nchar(label) > length_limit) {
        words <- strsplit(label, " ")[[1]]
        result <- ""
        current_line <- ""
        lines_used <- 1
        
        for (word in words) {
          test_line <- if (nchar(current_line) == 0) word else paste(current_line, word)
          if (nchar(test_line) > length_limit && lines_used < line_limit) {
            result <- if (nchar(result) == 0) current_line else paste(result, current_line, sep = "\n")
            current_line <- word
            lines_used <- lines_used + 1
          } else {
            current_line <- test_line
          }
        }
        
        # Add the last line
        if (nchar(current_line) > 0) {
          result <- if (nchar(result) == 0) current_line else paste(result, current_line, sep = "\n")
        }
        
        return(result)
      } else {
        return(label)
      }
    }, USE.NAMES = FALSE)
  }
  
  # Apply label wrapping
  dataset$term_name_wrapped <- wrap_labels(dataset$term_name)
  
  # Create plot
  ggplot(dataset) +
    geom_point(aes(x = (intersection_size/term_size)/(term_size/number_of_expressed_genes), 
                   y = term_name_wrapped, 
                   color = p_value)) + 
    labs(y = paste0("Top 20 ", GO_type, " terms"), 
         x = x_lab,
         color = "p.adjust") +
    ggtitle(dataset_name) +
    theme(axis.text.y = element_text(size = rel(0.8), hjust = 1),
          plot.title = element_text(hjust = 0.5))
}

plot_top_GO_refined = function(data, GO_type, dataset_name, 
                                 x_lab = "Gene Ratio", 
                                 min_point_size = 2,
                                 max_point_size = 8) {
  
  # Input validation
  if (!dataset_name %in% names(data)) stop("Dataset not found")
  dataset <- data[[dataset_name]] |> subset(source == GO_type)

  # Prepare data
  dataset <- dataset |> 
    mutate(
      gene_ratio = (intersection_size / term_size) / 
        (term_size / number_of_expressed_genes)
    ) |> 
    arrange(p_value) |> 
    head(20)
  
  # Determine reasonable breaks for gene count legend
  size_breaks <- pretty(dataset$intersection_size, n = 4)
  
  # Create plot
  ggplot(dataset, aes(x = gene_ratio, 
                      y = reorder(term_name, -p_value))) +
    geom_point(aes(size = intersection_size, color = p_value), alpha = 0.7) +
    scale_size_continuous(
      range = c(min_point_size, max_point_size),
      breaks = size_breaks,
      name = "Gene count"
    ) +
    scale_color_continuous(
      low = "red", high = "blue",
      trans = "log10",
      name = "p.adjust",
      labels = scales::scientific
    ) +
    labs(
      y = paste0("Top 20 ", GO_type, " terms"), 
      x = x_lab,
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_line(linetype = "dotted"),
      legend.position = "right",
      legend.box = "vertical"
    ) +
    guides(
      size = guide_legend(order = 1),
      color = guide_colorbar(order = 2)
    )
}



# Calculate semantic similarity matrices for each gene
top_GO_sim_matrix = function(data, GO_label, dataset_name){
  dataset <- subset(data[[dataset_name]], source == paste0('GO:',GO_label))
  simMatrix <- calculateSimMatrix(dataset$term_id,
                                  orgdb="org.Dm.eg.db",
                                  ont=GO_label,
                                  method="Rel")
}

# Reduced GO set based on semantic similarity matrices, tailored for sets where there are >20 GO terms
top_GO_reduced = function(simMatrix, GO_label, similarity, simMatrix_name){
reducedTerms <- reduceSimMatrix(simMatrix[[simMatrix_name]],
                                scores="uniqueness",
                                threshold=similarity,
                                orgdb="org.Dm.eg.db")
}

# GO three sets of genes (e.g., Ctrl unique, HS unique or both) and get the GO terms and gene overlap for each set
GO_Threeset_Venn = function(GO_data, gene_data, quantile, is_data_reduced, Var_or_Mean, GO_category){
  title <- textGrob(paste0(Var_or_Mean," Quantile ",quantile,' Genes - GO:',GO_category," Terms"),
                    gp = gpar(fontsize = 10)
  )
  Ctrlset <- paste0(Var_or_Mean,"CtrlQuantile",quantile)
  HSset <- paste0(Var_or_Mean,"HSQuantile",quantile)
  Bothset <- paste0(Var_or_Mean,"BothQuantile",quantile)
  # Get gene names
  Ctrl <- c(gene_data[[Ctrlset]],gene_data[[Bothset]])
  HS <- c(gene_data[[HSset]],gene_data[[Bothset]])
  # Get GO term names
  if (is_data_reduced==T){
    Ctrl_term <- GO_data[[Ctrlset]]$term
    HS_term <- GO_data[[HSset]]$term
    Both_term <- GO_data[[Bothset]]$term
  }
  else {
    Ctrl_term <- GO_data[[Ctrlset]]$term_name
    HS_term <- GO_data[[HSset]]$term_name
    Both_term <- GO_data[[Bothset]]$term_name
  }
  Venn <- ggvenn(list(Ctrl=Ctrl,HS=HS),
                 stroke_size = 0.5,
                 set_name_size = 5)
  text_A <- paste("Unique to Ctrl:","\n", paste(Ctrl_term, collapse = "\n"))
  text_B <- paste("Unique to HS:","\n", paste(HS_term, collapse = "\n"))
  text_AB <- paste("Both:","\n", paste(Both_term, collapse = "\n"))
  text_grobA <- textGrob(text_A,gp = gpar(fontsize = 7))
  text_grobB <- textGrob(text_B,gp = gpar(fontsize = 7))
  text_grobAB <- textGrob(text_AB,gp = gpar(fontsize = 7))
  return(grid.arrange(
    grobs = list(title, Venn, text_grobA,text_grobAB,text_grobB),
    heights = c(1,3,7),
    widths = c(1,1,1),
    layout_matrix = rbind(
      1,  # Title spans the top row
      2,  # Venn diagram spans the middle row
      c(3, 4, 5)   # Bottom elements in the third row
    ))
  )
}

#
# Function to convert list-type columns to strings for saving the full GO enrichment tables
sanitize_for_csv <- function(df) {
  list_cols <- sapply(df, is.list)
  if (any(list_cols)) {
    df[list_cols] <- lapply(df[list_cols], function(col) sapply(col, toString))
  }
  return(df)
}


# Function to convert into a clusterprofiler object 
# Allows a better visualisation (treeplot) to be carried out
gprofiler_to_enrichResult <- function(gost_result,
                                      organism = "UNKNOWN",
                                      keytype = "UNKNOWN",
                                      ontology = "UNKNOWN",
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.2,
                                      pAdjustMethod = "BH") {
  if (!"result" %in% names(gost_result)) {
    stop("Invalid g:Profiler result: must contain 'result' data frame.")
  }
  
  df <- gost_result$result
  
  if (nrow(df) == 0) {
    stop("g:Profiler result contains no enrichment terms.")
  }
  
  # Construct geneID column (genes associated with each term)
  df$geneID <- sapply(df$intersection, function(x) paste(x, collapse = "/"))
  
  # Create GeneRatio and BgRatio columns
  df$GeneRatio <- paste0(df$significant, "/", df$query_size)
  df$BgRatio <- paste0(df$effective_domain_size, "/", df$query_size)
  
  # Compose the result data frame
  enr_df <- df %>%
    dplyr::transmute(
      ID = term_id,
      Description = term_name,
      GeneRatio,
      BgRatio,
      pvalue = p_value,
      p.adjust = p_adj,
      qvalue = p_adj,  # or use df$p_value if no separate q-values
      geneID,
      Count = significant
    )
  
  # Flatten all input genes (used in constructing enrichResult)
  all_genes <- unique(unlist(df$intersection))
  
  # Create enrichResult object
  enrich_obj <- methods::new("enrichResult",
                             result = enr_df,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             qvalueCutoff = qvalueCutoff,
                             gene = all_genes,
                             universe = NULL,
                             geneSets = NULL,
                             organism = organism,
                             keytype = keytype,
                             ontology = ontology,
                             readable = FALSE)
  
  return(enrich_obj)
}
