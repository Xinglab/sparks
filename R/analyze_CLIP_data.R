###### ENRICHMENT TEST #####

#' Perform CLIP Enrichment Test for All regions individually
#'
#' This function performs test for enrichment in CLIP peaks in each region
#' relative to splicing invidiaully.
#' The function includes "event_anno" to define AS event junction type, whether
#' known AS events or novel combination of junctions or novel splice sites.
#' Also this function supports background setting via all_events.
#'
#'
#' @param object SPARKS object
#' @param as_events AS events of interest
#' @param event_type either SE, A3SS, A5SS
#' @param event_anno Junction information - KnownJC, NovelJC, NovelSS - option to limit
#' @param all_events Background + Foreground AS events to limit the scope of the analysis
#' @return Result df in 5 columns - "pval", "odds_ratio", "peak_count", "RBP", "Region"
#' @export
perform_CLIP_enrichment_test_regional <- function(object,
                                                  as_events,
                                                  event_type = "SE",
                                                  event_anno = c(),
                                                  all_events = c()){


  asso_df_region_list <- object@clip_result_list[[event_type]]
  if (event_type == "SE"){
    region_level <- c("Upstream_5ss_exon",
                      "Upstream_5ss_intron",
                      "Cassette_3ss_intron",
                      "Cassette_3ss_exon",
                      "Cassette_5ss_exon",
                      "Cassette_5ss_intron",
                      "Dnstream_3ss_intron",
                      "Dnstream_3ss_exon")
  } # TODO - implement parts for other splice types

  region_test_df <- lapply(region_level, function(region_of_interest){
    print(sprintf("Processing %s", region_of_interest))
    clip_peak_df <- asso_df_region_list[[region_of_interest]]
    if (length(event_anno) > 0){  # if AS junction types are selected
      clip_peak_df <- clip_peak_df[object@exon_annotation$SE[rownames(clip_peak_df), , drop=FALSE]$annotation %in% event_anno, ]
    }
    if (length(all_events) > 0){  # if background are selected
      clip_peak_df <- clip_peak_df[all_events, ]
    }

    test_result_new <- lapply(colnames(clip_peak_df), function(test_col){
      fisher_test_data <- data.frame(peak = clip_peak_df[[test_col]],
                                     asso = ifelse(rownames(clip_peak_df) %in% as_events, "Asso", "No_asso"),
                                     row.names = rownames(clip_peak_df))

      tryCatch(test_result <<- fisher.test(table(fisher_test_data), alternative = "less"), error = function(e){test_result <<- NULL})
      if (is.null(test_result)){  # put in default value if test result is null
        p_val <- 1
        odds_ratio <- 1
      } else {
        p_val <- test_result$p.value
        odds_ratio <- test_result$estimate
      }

      if (dim(table(fisher_test_data))[1] == c(2, 2)){  # if it is well-formed contingency table
        ratio_str <- sprintf("%s/%s", table(fisher_test_data)[2, 1], sum(table(fisher_test_data)[, 1]))
      } else {
        ratio_str <- "NA"
      }

      # append the result
      return(c(p_val, odds_ratio, ratio_str))
    })
    test_df <- as.data.frame(do.call(rbind, test_result_new))
    colnames(test_df) <- c("pval", "odds_ratio", "peak_count")
    test_df$study_rbp <- colnames(clip_peak_df)
    test_df$region <- region_of_interest
    return(test_df)
  })

  region_pval_df <- do.call(rbind, region_test_df)
  region_pval_df$pval <- as.numeric(region_pval_df$pval)
  return(region_pval_df)
}


#' Perform CLIP Enrichment Test for All regions individually for directionality
#'
#' This function performs test for enrichment in CLIP peaks in each region
#' relative to splicing invidiaully.
#' This function is designed to test directional enrichment, given positive and negative events
#'
#'
#'
#' @param object SPARKS object
#' @param pos_events AS events of interest in one group
#' @param event_anno AS events of interest in the other group
#' @param event_type either SE, A3SS, A5SS
#' @return Result df in 5 columns - "pval", "odds_ratio", "peak_count", "RBP", "Region"
#' @export
perform_CLIP_enrichment_test_direction <- function(object, pos_events, neg_events, event_type = "SE"){
  # define background
  events_of_interest <- unique(c(pos_events, neg_events))
  # extract relevant clip result
  asso_df_region_list <- object@clip_result_list[[event_type]]
  if (event_type == "SE"){
    region_level <- c("Upstream_5ss_exon",
                      "Upstream_5ss_intron",
                      "Cassette_3ss_intron",
                      "Cassette_3ss_exon",
                      "Cassette_5ss_exon",
                      "Cassette_5ss_intron",
                      "Dnstream_3ss_intron",
                      "Dnstream_3ss_exon")
  } # TODO - implement parts for other splice types

  # run the analysis seperately
  region_test_df <- lapply(region_level, function(region_of_interest){
    print(sprintf("Processing %s", region_of_interest))
    clip_peak_df <- asso_df_region_list[[region_of_interest]]

    # take subset
    clip_peak_df <- clip_peak_df[events_of_interest, ]

    test_result_new <- lapply(colnames(clip_peak_df), function(test_col){
      fisher_test_data <- data.frame(peak = clip_peak_df[[test_col]],
                                     asso = ifelse(rownames(clip_peak_df) %in% pos_events, "Positive Asso", "Negative Asso"),
                                     row.names = rownames(clip_peak_df))

      tryCatch(test_result <- fisher.test(table(fisher_test_data)), error = function(e){test_result <<- "Not_Found"})
      if (test_result == "Not_Found"){
        p_val <- 1
        odds_ratio <- 1
      } else {
        p_val <- test_result$p.value
        odds_ratio <- test_result$estimate
      }

      if (dim(table(fisher_test_data))[1] == c(2,2)){ # if it is well-formed contingency table
        ratio_str <- sprintf("%s/%s vs. %s/%s",
                             table(fisher_test_data)['peak', 'Positive Asso'],
                             sum(table(fisher_test_data)[,'Positive Asso']),
                             table(fisher_test_data)['peak', 'Negative Asso'],
                             sum(table(fisher_test_data)[,'Negative Asso'])
        )
      } else {
        ratio_str <- "NA"
      }

      # append the result
      return(c(p_val, odds_ratio, ratio_str))
    })
    test_df <- as.data.frame(do.call(rbind, test_result_new))
    colnames(test_df) <- c("pval", "odds_ratio", "peak_count")
    test_df$study_rbp <- colnames(clip_peak_df)
    test_df$region <- region_of_interest
    test_df$pval <- as.numeric(test_df$pval)

    return(test_df)
  })



  region_pval_df <- do.call(rbind, region_test_df)
  return(region_pval_df)
}




#' Perform CLIP Enrichment Test for All regions merged
#'
#' This function performs test for enrichment in CLIP peaks in all regions combined.
#' Unlike the other function, it considers the peak if it is in any of the regions in splicing area.
#' The function includes "event_anno" to define AS event junction type, whether
#' known AS events or novel combination of junctions or novel splice sites.
#' Also this function supports background setting via all_events.
#'
#'
#' @param object SPARKS object
#' @param as_events AS events of interest
#' @param event_type either SE, A3SS, A5SS
#' @param event_anno Junction information - KnownJC, NovelJC, NovelSS - option to limit
#' @param all_events Background + Foreground AS events to limit the scope of the analysis
#' @return Result df in 5 columns - "pval", "odds_ratio", "peak_count", "RBP", "Region"
#' @export
perform_CLIP_enrichment_test_any_region <- function(object, as_events, event_type = "SE", event_anno = c(), all_events = c()){


  asso_df_region_list <- object@clip_result_list[[event_type]]
  if (event_type == "SE"){
    region_level <- c("Upstream_5ss_exon",
                      "Upstream_5ss_intron",
                      "Cassette_3ss_intron",
                      "Cassette_3ss_exon",
                      "Cassette_5ss_exon",
                      "Cassette_5ss_intron",
                      "Dnstream_3ss_intron",
                      "Dnstream_3ss_exon")
  } # TODO - implement parts for other splice types
  clip_3d_matrix <- abind(asso_df_region_list, along = 3)
  clip_experiment_list <- colnames(clip_3d_matrix[, ,1])

  # calculate whether entire AS events have any CLIP peaks
  sum_clip_df <- do.call(cbind, lapply(clip_experiment_list, function(clip_experiment){
    clip_matrix <- clip_3d_matrix[, clip_experiment, ]

    clip_peak_count <- rowSums(clip_matrix == "peak", na.rm = T)

    clip_status <- ifelse(clip_peak_count > 0, "peak", "no_peak")
    return(clip_status)
  }))
  colnames(sum_clip_df) <- clip_experiment_list

  clip_peak_df <- sum_clip_df
  if (length(event_anno) > 0){
    clip_peak_df <- clip_peak_df[object@exon_annotation$SE[rownames(clip_peak_df), , drop=FALSE]$annotation %in% event_anno, ]
  }
  if (length(all_events) > 0){
    clip_peak_df <- clip_peak_df[all_events, ]
  }

  test_result_new <- lapply(colnames(clip_peak_df), function(test_col){
    fisher_test_data <- data.frame(peak = clip_peak_df[,test_col],
                                   asso = ifelse(rownames(clip_peak_df) %in% as_events, "Asso", "No_asso"),
                                   row.names = rownames(clip_peak_df))

    tryCatch(test_result <<- fisher.test(table(fisher_test_data), alternative = "less"), error = function(e){test_result <<- NULL})
    if (is.null(test_result)){
      p_val <- 1
      odds_ratio <- 1
    } else {
      p_val <- test_result$p.value
      odds_ratio <- test_result$estimate
    }

    if (dim(table(fisher_test_data))[1] == c(2,2)){ # if it is well-formed contingency table
      ratio_str <- sprintf("%s/%s", table(fisher_test_data)[2,1], sum(table(fisher_test_data)[,1]))
    } else {
      ratio_str <- "NA"
    }

    # append the result
    return(c(p_val, odds_ratio, ratio_str))
  })
  test_df <- as.data.frame(do.call(rbind, test_result_new))
  colnames(test_df) <- c("pval", "odds_ratio", "peak_count")
  test_df$study_rbp <- colnames(clip_peak_df)
  test_df$region <- "combined"

  region_pval_df <- test_df
  region_pval_df$pval <- as.numeric(region_pval_df$pval)

  return(region_pval_df)
}




#' Perform CLIP Enrichment Test for All regions merged in directinoal
#'
#' This function performs test for enrichment in CLIP peaks in all regions combined.
#' Unlike the other function, it considers the peak if it is in any of the regions in splicing area.
#' The function includes "event_anno" to define AS event junction type, whether
#' known AS events or novel combination of junctions or novel splice sites.
#' Also this function supports background setting via all_events.
#'
#'
#' @param object SPARKS object
#' @param as_events AS events of interest
#' @param event_type either SE, A3SS, A5SS
#' @param event_anno Junction information - KnownJC, NovelJC, NovelSS - option to limit
#' @param all_events Background + Foreground AS events to limit the scope of the analysis
#' @return Result df in 5 columns - "pval", "odds_ratio", "peak_count", "RBP", "Region"
#' @export
perform_CLIP_enrichment_test_any_region_directional <- function(object, pos_events, neg_events, event_type = "SE"){


  asso_df_region_list <- object@clip_result_list[[event_type]]
  if (event_type == "SE"){
    region_level <- c("Upstream_5ss_exon",
                      "Upstream_5ss_intron",
                      "Cassette_3ss_intron",
                      "Cassette_3ss_exon",
                      "Cassette_5ss_exon",
                      "Cassette_5ss_intron",
                      "Dnstream_3ss_intron",
                      "Dnstream_3ss_exon")
  } # TODO - implement parts for other splice types
  clip_3d_matrix <- abind(asso_df_region_list, along = 3)
  clip_experiment_list <- colnames(clip_3d_matrix[, ,1])

  # calculate whether entire AS events have any CLIP peaks
  sum_clip_df <- do.call(cbind, lapply(clip_experiment_list, function(clip_experiment){
    clip_matrix <- clip_3d_matrix[, clip_experiment, ]

    clip_peak_count <- rowSums(clip_matrix == "peak", na.rm = T)

    clip_status <- ifelse(clip_peak_count > 0, "peak", "no_peak")
    return(clip_status)
  }))
  colnames(sum_clip_df) <- clip_experiment_list

  # subset the events
  clip_peak_df <- sum_clip_df[rownames(sum_clip_df) %in% c(pos_events, neg_events), ]

  test_result_new <- lapply(colnames(clip_peak_df), function(test_col){
    fisher_test_data <- data.frame(peak = clip_peak_df[, test_col],
                                   asso = ifelse(rownames(clip_peak_df) %in% pos_events, "Positive Asso", "Negative Asso"),
                                   row.names = rownames(clip_peak_df))

    tryCatch(test_result <- fisher.test(table(fisher_test_data)), error = function(e){test_result <<- "Not_Found"})
    if (test_result == "Not_Found"){
      p_val <- 1
      odds_ratio <- 1
    } else {
      p_val <- test_result$p.value
      odds_ratio <- test_result$estimate
    }

    if (dim(table(fisher_test_data))[1] == c(2,2)){ # if it is well-formed contingency table
      ratio_str <- sprintf("%s/%s vs. %s/%s",
                           table(fisher_test_data)['peak', 'Positive Asso'],
                           sum(table(fisher_test_data)[,'Positive Asso']),
                           table(fisher_test_data)['peak', 'Negative Asso'],
                           sum(table(fisher_test_data)[,'Negative Asso'])
      )
    } else {
      ratio_str <- "NA"
    }

    # append the result
    return(c(p_val, odds_ratio, ratio_str))
  })
  test_df <- as.data.frame(do.call(rbind, test_result_new))
  colnames(test_df) <- c("pval", "odds_ratio", "peak_count")
  test_df$study_rbp <- colnames(clip_peak_df)
  test_df$region <- "combined"

  region_pval_df <- test_df
  region_pval_df$pval <- as.numeric(region_pval_df$pval)

  return(region_pval_df)
}



##### QUERY FOR OTHER ANALYSIS #####

#' Query CLIP peak df for given RBP
#'
#'
#' @param object SPARKS object
#' @param target Specific RBPs to query CLIP peaks for
#' @return CLIP peak dataframe - columns are region, rows are AS events, elements are whether it has peak at that position or not
#' @export
query_clip_peak_df_for_RBP <- function(object, target, event_type = "SE"){
  # call peak df for specific RBP
  peak_info_df <- do.call(cbind, lapply(names(object@clip_result_list[[event_type]]), function(x){
    target_peak_info <- object@clip_result_list[[event_type]][[x]][, target, drop = FALSE]
    # target_peak_info <- ica_object@clip_result_list$SE[[x]][aaa_as_events, target, drop = FALSE]

    colnames(target_peak_info) <- x
    return(target_peak_info)
  }))
  return(peak_info_df)
}


#' Query CLIP peak df for given AS events
#'
#'
#' @param object SPARKS object
#' @param as_events AS events to get query for
#' @return CLIP peak dataframe - columns are region, rows are AS events, elements are whether it has peak at that position or not
#' @export
query_clip_peak_df_for_AS_events <- function(object, as_events, event_type = "SE"){
  asso_df_region_list <- object@clip_result_list[[event_type]]
  if (event_type == "SE"){
    region_level <- c("Upstream_5ss_exon",
                      "Upstream_5ss_intron",
                      "Cassette_3ss_intron",
                      "Cassette_3ss_exon",
                      "Cassette_5ss_exon",
                      "Cassette_5ss_intron",
                      "Dnstream_3ss_intron",
                      "Dnstream_3ss_exon")
  } # TODO - implement parts for other splice types
  clip_3d_matrix <- abind(asso_df_region_list, along = 3)
  clip_experiment_list <- colnames(clip_3d_matrix[, ,1])

  # calculate whether entire AS events have any CLIP peaks
  sum_clip_df <- do.call(cbind, lapply(clip_experiment_list, function(clip_experiment){
    clip_matrix <- clip_3d_matrix[, clip_experiment, ]

    clip_peak_count <- rowSums(clip_matrix == "peak", na.rm = T)

    clip_status <- ifelse(clip_peak_count > 0, "peak", "no_peak")
    return(clip_status)
  }))
  colnames(sum_clip_df) <- clip_experiment_list

  # subset the events of interest
  sum_clip_df_eoi <- as.data.frame(sum_clip_df[as_events[as_events %in% rownames(sum_clip_df)], ])

  return(sum_clip_df_eoi)
}




#' Generate CLIP Peak enrichment HEATMAP
#'
#'
#' @param region_pval_select result df from CLIP peak enrichment test
#' @param event_type either SE, A3SS, A5SS
#' @return heatmap generated by ggplot
#' @export
generate_CLIP_enrichment_heatmap <- function(region_pval_select,
                                             event_type = "SE"){
  if (event_type == "SE"){
    region_level <- c("Upstream_5ss_exon",
                      "Upstream_5ss_intron",
                      "Cassette_3ss_intron",
                      "Cassette_3ss_exon",
                      "Cassette_5ss_exon",
                      "Cassette_5ss_intron",
                      "Dnstream_3ss_intron",
                      "Dnstream_3ss_exon")
  }
  # re-order x-axis (regions relative to exon)
  region_pval_select$region <- factor(region_pval_select$region, levels = region_level)
  region_pval_select$study_rbp <- factor(region_pval_select$study_rbp,
                                         levels =  rev(unique(region_pval_select$study_rbp[order(unlist(lapply(region_pval_select$study_rbp,
                                                                                                               function(x) strsplit(x, "_")[[1]][2])))])))
  p <- ggplot(region_pval_select, aes(x = region, y = study_rbp, fill = -log10(as.numeric(pval)))) +
    # geom_raster(aes(alpha = -log2(as.numeric(odds_ratio)))) +
    geom_raster() +
    geom_text(aes(label = peak_count)) +
    scale_fill_distiller(palette = "Reds", direction = 1,
                         limits = c(0, max(-log10(as.numeric(region_pval_select$pval))))) +
    theme(legend.position = "bottom") +
    scale_x_discrete(expand = c(0, 0))
  return(p)
}



#' Generate CLIP peak annotated DS events
#'
#'
#' @param object SPARKS object
#' @param target Target RBP
#' @param event_type either SE, A3SS, A5SS
#' @return result dataframe with clip peak presence and beta from MATS
#' @export
call_CLIP_targets_for_ENCODE_data <- function(object, target, event_type = "SE"){


  # call known JC events
  known_jc <- rownames(subset(object@exon_annotation[[event_type]], annotation == "Known_JC"))

  # call DS events
  mats_df <- object@MATS_list[[event_type]]
  mats_df_known <- mats_df[rownames(mats_df) %in% known_jc, ]

  mats_df_sig <- subset(mats_df_known, abs(beta) > 0.1 & fdr < 0.1)
  ds_events <- rownames(mats_df_sig)

  # check if the target has corresponding CLIP peaks
  if (target %in% colnames(object@clip_result_list$SE$Cassette_3ss_intron)){
    # call CLIP peaks for corresponding RBPs
    clip_df <- query_clip_peak_df_for_RBP(object, target, event_type = event_type)

    # flatten by presence of peak in any region
    clip_df_sum <- rowSums(clip_df == "peak")

    clip_sig_events <- names(clip_df_sum[clip_df_sum > 0])


    # generate combined plot df
    result_df <- data.frame(event = ds_events,
                            target = ifelse(ds_events %in% clip_sig_events, "primary", "secondary"),
                            beta = mats_df_sig$beta,
                            fdr = mats_df_sig$fdr)
  } else {
    result_df <- data.frame(event = ds_events,
                            target = "Not_determined",
                            beta = mats_df_sig$beta,
                            fdr = mats_df_sig$fdr)
  }
  return(result_df)
}

