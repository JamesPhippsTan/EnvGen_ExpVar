make_quantiles <-function(data,quantile_column,n_quantiles,name_of_quantiles) {
  quantile <- ceiling(nrow(data)/n_quantiles)
  data <- data[order(quantile_column),]
  for (i in 1:n_quantiles){
    firstrow <- 1 + quantile*(i-1)
    if (i < n_quantiles){
      lastrow <- quantile*i
    }
    else {
      lastrow <- nrow(data)
    }
    data[rownames(data)[firstrow:lastrow],name_of_quantiles] <- i
  }
  data[[name_of_quantiles]] <- as.factor(data[[name_of_quantiles]])
  return(data)
}

make_rankplot_df <- function(df_with_ranks){
  rankplot_df <- data.frame( 
    Gene = c(rownames(df_with_ranks),rownames(df_with_ranks)),
    Diet = rep(c("Ctrl", "HS"), each = nrow(df_with_ranks)),
    Variability_Rank = c(df_with_ranks[["Ctrl_Variability_Rank"]],  df_with_ranks[["HS_Variability_Rank"]]),
    Mean_Rank = c(df_with_ranks[["Ctrl_Mean_Rank"]],  df_with_ranks[["HS_Mean_Rank"]])
  )
  return(rankplot_df)
}



percent_per_quantile <-function(data,column_with_quantiles,column_to_percent_name,new_column_name){
  column_with_quantiles <- as.numeric(column_with_quantiles)
  quantiles  <- sort(unique(column_with_quantiles))
  last_quantile <- max(quantiles)
  genes_per_quantile <- ceiling(nrow(data)/last_quantile)
  data <- data[order(column_with_quantiles),]
  for (i in quantiles){
    firstrow <- 1 + genes_per_quantile*(i-1)
    if (i < last_quantile){
      lastrow <- genes_per_quantile*i
    }
    else {
      lastrow <- nrow(data)
    }
    quantile_subset <- subset(data,column_with_quantiles == i)
    data[rownames(data)[firstrow:lastrow],new_column_name] <- sum(quantile_subset[,column_to_percent_name])/genes_per_quantile
  }
  return(data)
}
