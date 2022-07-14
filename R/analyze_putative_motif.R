

##### SUBFUNCTIONS #####

cleanup_data_for_fisher <- function(peak_df, region, mut_asso_events, spl_event_type){
  # binarize binding site count
  peak_count <- peak_df
  peak_count$region_ind <- ifelse(as.numeric(peak_count$peak_count) > 0, "RBP BS", "No RBP BS")
  # append splicing type
  # TODO - this is fixed. fix it for later use
  # peak_count$id_fixed <- unlist(lapply(rownames(peak_count), function(x) sprintf("%s:%s", x, spl_event_type)))

  # set up positive and negative group separately
  peak_count$asso_ind <- ifelse(rownames(peak_count) %in% mut_asso_events, "Associated with Mutation", "Not associated with Mutation")
  # print(peak_count)
  # prepare data for fisher test
  annotated_count_df <- peak_count[, c("asso_ind", "region_ind")]
  return(annotated_count_df)
}

cleanup_data_for_fisher_directional <- function(peak_df, region, pos_asso_events, spl_event_type){
  # binarize binding site count
  peak_count <- peak_df
  peak_count$region_ind <- ifelse(as.numeric(peak_count$peak_count) > 0, "RBP BS", "No RBP BS")
  # append splicing type
  # TODO - this is fixed. fix it for later use
  # peak_count$id_fixed <- unlist(lapply(rownames(peak_count), function(x) sprintf("%s:%s", x, spl_event_type)))

  # set up positive and negative group separately
  peak_count$asso_ind <- ifelse(rownames(peak_count) %in% pos_asso_events, "Positively Associated", "Negatively Associated")
  # print(peak_count)
  # prepare data for fisher test
  annotated_count_df <- peak_count[, c("asso_ind", "region_ind")]
  return(annotated_count_df)
}


calculate_fisher_exact_stat <- function(peak_df, region, mut_asso_events, spl_event_type){

  annotated_count_df <- cleanup_data_for_fisher(peak_df, region, mut_asso_events, spl_event_type)
  # perform test
  tryCatch(test_result <- fisher.test(table(annotated_count_df), alternative = "less"), error = function(e){test_result <<- NULL})
  if (is.null(test_result)){
    p_val <<- 1
    odds_ratio <<- 1
  } else {
    p_val <<- test_result$p.value
    odds_ratio <<- test_result$estimate
  }

  if (dim(table(annotated_count_df))[1] == 2){ # if it is well-formed contingency table
    print(table(annotated_count_df))
    ratio_str <- sprintf("%s/%s", table(annotated_count_df)[1,2], sum(table(annotated_count_df)[1,]))
  } else {
    ratio_str <- "NA"
  }

  # append the result
  return(c(p_val, odds_ratio, ratio_str))
}


calculate_fisher_exact_stat_directional <- function(peak_df, region, mut_asso_events, spl_event_type){

  annotated_count_df <- cleanup_data_for_fisher_directional(peak_df, region, mut_asso_events, spl_event_type)
  # perform test
  tryCatch(test_result <- fisher.test(table(annotated_count_df)), error = function(e){test_result <<- NULL})
  if (is.null(test_result)){
    p_val <<- 1
    odds_ratio <<- 1
  } else {
    p_val <<- test_result$p.value
    odds_ratio <<- test_result$estimate
  }

  if (dim(table(annotated_count_df))[1] == 2 & dim(table(annotated_count_df))[2] == 2){ # if it is well-formed contingency table
    print(table(annotated_count_df))
    ratio_str <- sprintf("%s/%s vs. %s/%s",
                         table(annotated_count_df)['Positively Associated', 'RBP BS'],
                         sum(table(annotated_count_df)['Positively Associated', ]),
                         table(annotated_count_df)['Negatively Associated', 'RBP BS'],
                         sum(table(annotated_count_df)['Negatively Associated', ]))
  } else {
    ratio_str <- "NA"
  }

  # append the result
  return(c(p_val, odds_ratio, ratio_str))
}


##### MAIN FUNCTIONS #####
perform_enrichment_analysis <- function(object, mutation_event, spl_event_type = "SE"){
  print(sprintf("Performing Motif Enrichment Analysis for %s", spl_event_type))
  # select the matrix list for event type
  count_matrix_list <- object@motif_enrichment_df[[spl_event_type]]

  # extract associated splicing event
  mutation_asso_event_df <- object@linear_model_pval[[spl_event_type]][mutation_event, ]
  mutation_asso_event <- colnames(mutation_asso_event_df)[mutation_asso_event_df > 0]
  regional_result <- list()
  for (region in names(count_matrix_list)){
    print(sprintf("- performing analysis for %s region", region))
    count_matrix <- count_matrix_list[[region]]

    result_list <- list()
    dummy <- apply(count_matrix, 2, function(peak_count){
      peak_df <- as.data.frame(peak_count)
      rownames(peak_df) <- rownames(count_matrix)
      # print(peak_df)
      test_result_a <- calculate_fisher_exact_stat(peak_df, region, mutation_asso_event, spl_event_type)
      # print(test_result)
      neglogpval <- -log10(as.numeric(test_result_a[[1]][1]))
      odds_ratio <- -log2(as.numeric(test_result_a[[2]][1]))
      ratio_str <- test_result_a[[3]][1]
      result_df <- data.frame("Region" = region,
                              "neglogpval" = neglogpval,
                              "Odds_ratio" = odds_ratio,
                              "Motif_Occurrence" = ratio_str)
      result_list[[length(result_list) + 1]] <<- result_df
    })
    # gather and rename the results
    result_df <- do.call(rbind, result_list)
    result_df$RBP <- colnames(count_matrix)

    # add this to outside list to gather results for all three regions
    regional_result[[length(regional_result) + 1]] <- result_df
  }
  analayis_df <- do.call(rbind, regional_result)

  # store the data
  object@motif_analysis_result[[mutation_event]] <- analayis_df
  return(object)

}



