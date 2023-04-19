##### SUBFUNCTIONS #####


#' Calculate adjusted p-value for SPARKS analysis results
#'
#' This function sets the padj.
#' Any signatures with the same sign between positive ES and nevative ES are set
#' to have padj of 1, so they wouldn't be considered significant.
#' This approach is similar to C-MAP paper (Lamb et al., Science 2006) that
#' they remove any signatures with same sign are removed.
#'
#' @param sparks_result SPARKS analysis dataframe
#' @param sig_test_method p.adjust method - "BH" for FDR, "bonferroni" (default) for FWER
# #' @param same_sign_ES_threshold Set any signatures with abs score below this when sign from pos and neg signature sets are the same
#'
#' @return SPARKS analysis dataframe with padj added
#' @export
calculate_SPARKS_padj <- function(sparks_result, sig_test_method = "bonferroni", same_sign_ES_threshold = 0.5){
  # calculate raw padj
  sparks_result$padj <- p.adjust(sparks_result$gsea_combined_pval,
                                 method = sig_test_method)

  # flatten if the sign is the same but the score is below the threshold
  # - replacement needs at least one contradicting entry, otherwise it throws error
  # - So, we check before we do this
  if (sum(sign(sparks_result$gsea_pos_score) == sign(sparks_result$gsea_neg_score)) > 0){
    sparks_result[(sign(sparks_result$gsea_pos_score) == sign(sparks_result$gsea_neg_score)) &
                    (sparks_result$gsea_score_abs < same_sign_ES_threshold), ]$padj <- 1

  }
  return(sparks_result)
}


#' Calculate enrichment score using FGSEA
#'
#' @param gsea_library Signature event list (positive and negative) from `extract_GSEA_significant_events` function
#' @param study_rank AS events sorted by delta PSI
#'
#' @return gsea_result_df Dataframe with enrichment score and pval for positive and negative event sets
#' @export
calculate_fgsea_score <- function(gsea_library, study_rank){
  # calculate score using fgsea
  gsea_result <- fgsea::fgsea(pathways = gsea_library ,
                              stats = study_rank,
                              nproc = 1)

  # return 0 values if the method fails to run on either pos or neg
  if(dim(gsea_result)[1] < 2){
    output_df <- data.frame(score = 0,
                            pos_score = 0,
                            neg_score = 0,
                            pos_pval = 0,
                            neg_pval = 0,
                            pval = 1)
    return(output_df)
  }


  # gather the data

  # combine score
  neg_score <- subset(gsea_result, pathway == "negative")$ES
  pos_score <- subset(gsea_result, pathway == "positive")$ES

  combined_score <- pos_score - neg_score

  # combine pval
  neg_pval <- subset(gsea_result, pathway == "negative")$pval
  pos_pval <- subset(gsea_result, pathway == "positive")$pval
  # print(neg_pval)
  # print(gsea_result)

  if (is.na(neg_pval)){  # if we don't get the enrichment score for some reason..
    combined_pval <- metap::sumlog(c(pos_pval, 1))$p

  } else if (is.na(pos_pval)){
    combined_pval <- metap::sumlog(c(neg_pval, 1))$p

    # if(is.na(neg_pval) | is.na(pos_pval)){
    #   combined_pval <- 1
    # } else if (sign(pos_score) == sign(neg_score)){  # if the sign is the same, p-value should be 1
    #   # this is based on CMAP paper (Lamb et al., Science 2006)
    #   # the CMAP paper makes the score 0, but for our analysis
    #   # we handle this at the multiple testing correction step
  } else {
    fisher_result <- metap::sumlog(c(neg_pval, pos_pval))
    combined_pval <- fisher_result$p
  }

  output_df <- data.frame(score = combined_score,
                          pos_score = pos_score,
                          neg_score = neg_score,
                          pos_pval = pos_pval,
                          neg_pval = neg_pval,
                          pval = combined_pval)
  return(output_df)
}


#' Extract signature events for SPARKS
#'
#' This function extracts events with abs(delPSI) > 0.1 and store in a list, sepearated in
#' positive and negative set of AS events
#'
#' @param study_mats filtered MATS df
#'
#' @return event_list List of events with positive (delPSI > 0.1) and negative (delPSI < -0.1) AS changes
#' @export
extract_GSEA_significant_events <- function(study_mats){
  # divide pos events and neg events for query
  study_sig <- subset(study_mats, abs(beta) > 0.1)
  rownames(study_sig) <- study_sig$event

  study_beta_df <- study_sig[, c("beta"), drop = F ]
  study_beta_complete <- study_beta_df[!is.na(study_beta_df), , drop = F]
  colnames(study_beta_complete) <- 'beta'

  test_pos_events <- rownames(subset(study_beta_complete, beta > 0))
  test_neg_events <- rownames(subset(study_beta_complete, beta < 0))

  event_list <- list()
  event_list$positive <- test_pos_events
  event_list$negative <- test_neg_events
  return(event_list)
}


#' Rewrite event coordinates
#'
#' This function removes unnecessary info in AS event coordinates.
#'
#' @param x AS event
#'
#' @return filtered_AS_event
#' @export
rewrite_event_coordinates <- function(x){  # clean coordinates
  # split into elements for processing
  elements <- strsplit(x, ":")[[1]]
  # print(elements)

  # infer spl type
  spl_type <- elements[[11]]

  # process each splicing type info to extract only necessary info
  if (spl_type == "SE"){  # SE don't have directionality
    new_event <- paste0(c(elements[2:4], elements[8], elements[5:6], elements[9], elements[11]), collapse = ":")
  } else if (spl_type == "A3SS"){  # directionality matters for A3SS
    if (elements[[4]] == '+'){  # if positive
      new_event <- paste0(c(elements[2:4], elements[10], elements[5], elements[7], elements[11]), collapse = ":")
    } else if (elements[[4]] == "-"){  # if negative
      new_event <- paste0(c(elements[2:4], elements[9], elements[6], elements[8], elements[11]), collapse = ":")
    } else {
      stop()
    }
  } else if (spl_type == "A5SS"){
    if (elements[[4]] == '+'){  # if positive
      new_event <- paste0(c(elements[2:4], elements[8], elements[6], elements[9], elements[11]), collapse = ":")
    } else if (elements[[4]] == "-"){  # if negative
      new_event <- paste0(c(elements[2:4], elements[7], elements[5], elements[10], elements[11]), collapse = ":")
    } else {
      stop()
    }
  }
  return(new_event)
}


##### FUNCTION #####

#' Import MATS for SPARKS analysis from SPARKS object
#'
#' This function reads MATS file generated from the pipeline, keep known events,
#' and filter out events with the read counts below the count_threshold in any
#' samples.
#'
#' It also keeps only necessary information for the AS events defining junctions.
#'
#' @param input_sparks Input SPARKS object
#'
#' @param spl_type Optional splicing type in c("SE", "A3SS", "A5SS"), default = "SE"
#' @param count_threshold Average read count for AS events. Events with read counts below this number in average is removed
#' @param min_threshold Minimum read count for AS events. Events with read counts below this number in any sample is removed
#'
#' @return filtered_MATS_df
#' @export
import_SPARKS_MATS_for_analysis <- function(input_sparks,
                                            spl_type = "SE",
                                            count_threshold = 20,
                                            min_threshold = 20){
  # check if the event name is already processed
  if(length(strsplit(input_sparks@MATS_list[[spl_type]][1, 'event'], ":")[[1]]) == 11){  # if not processed
    # process the data for analysis
    known_events <- rownames(subset(input_sparks@exon_annotation[[spl_type]], annotation == "Known_JC"))
    study_mats_raw <- input_sparks@MATS_list[[spl_type]]

    # keep only the known AS events
    study_mats_known <- study_mats_raw[study_mats_raw$event %in% known_events, ]

    # update the event to short form
    study_mats_known$event <- unlist(lapply(study_mats_known$event,
                                            function(x) rewrite_event_coordinates(x)))
  } else {  # means the MATS is already processed, so no need to anything
    study_mats_known <- input_sparks@MATS_list[[spl_type]]
  }

  if ('pulled_delta_psi' %in% colnames(study_mats_known)){  # if old processed data with typo in pooled
    study_mats_known$beta <- as.numeric(study_mats_known$pulled_delta_psi)

    # change the column names form pulled to pooled
    colnames(study_mats_known) <- gsub("pulled", "pooled", colnames(study_mats_known))

  } else {  # if the column name is corrected
    study_mats_known$beta <- as.numeric(study_mats_known$pooled_delta_psi)
  }


  # calculate minimum_count
  study_mats_known$min_count <- unlist(lapply(study_mats_known$count_values,
                                              function(input_counts) min(do.call(as.numeric, strsplit(input_counts, ",")))))

  # filter by minimum count
  study_mats_temp <- subset(study_mats_known,
                            avg_count >= count_threshold &
                              min_count >= min_threshold) %>% arrange(-beta)

  return(study_mats_temp)
}





#' Run SPARKS analysis
#'
#' @param study_mats filtered MATS dataframe
#' @param kd_library Signature library
#' @param study Study name - this will be stored in "S2" column in the output
#' @param num_cores Number of cores for parallel run - default = 3
#' @param overlap_ratio_threshold Threshold to filter out AS events based on their occurrences in the library - default = 30%
#' @param library_list If given, SPARKS analysis will be performed in only this signatures
#' @param event_count_df If given, SPARKS will use this to enumerate events to filter out, rather than iterating the library
#' @param library_size Library size - currently this is 675. This is built in case for future cases where the signature library is expanded. Do not change this for now.
#'
#' @return SPARKS_analysis_df
#' @export
perform_SPARKS_analysis_with_overlap_filter <- function(study_mats,
                                                        kd_library,
                                                        study,
                                                        num_cores = 3,
                                                        overlap_ratio_threshold = 0.30,
                                                        library_list = c(),
                                                        event_count_df = data.frame(),
                                                        library_size = 675){
  method = "GSEA"

  ## remove frequently affected events by knockdown
  if (dim(event_count_df)[1] != 0){  # use pre-compiled list if it is there
    uninformative_events <- event_count_df[event_count_df$count > library_size * overlap_ratio_threshold, ]$event
  } else {  # calculate the count for each event and filter them
    #  events
    event_list <- list()

    # extract all sig events in the library
    dummy <- lapply(names(kd_library), function(experiment){
      sig_events <- unlist(extract_GSEA_significant_events(kd_library[[experiment]]))
      event_list[[length(event_list) + 1]] <<- sig_events
      return()
    })

    # generate list of AS events and count for them
    # null_event <- unique(unlist(event_list))
    event_count_df <- as.data.frame(table(unlist(event_list)), stringsAsFactors = FALSE)
    colnames(event_count_df) <- c("event", "count")

    # define uninformative events
    # uninformative_events <- names(null_weight)[null_weight > length(kd_library) * overlap_ratio_threshold]
    uninformative_events <- event_count_df[event_count_df$count > length(kd_library) * overlap_ratio_threshold, ]$event
  }

  # make dataframe ready for FGSEA
  study_mats_clean <- study_mats[, c("event", "pooled_delta_psi")]
  study_rank <- study_mats_clean$pooled_delta_psi
  names(study_rank) <- study_mats_clean$event

  # limit the scope of kd library if given
  if(length(library_list) > 0){
    kd_library <- kd_library[library_list]
  }

  # calculate correlation for signature
  test_result_df_list <- pbmcapply::pbmclapply(names(kd_library), function(signature){
    # filter out by count
    test_mats <- kd_library[[signature]]
    test_mats_filtered <- subset(test_mats, !(event %in% uninformative_events))

    result_df <- calculate_RBP_KD_correlation_from_mats(test_mats_filtered,  # apply this filter to here too
                                                        study_mats,
                                                        signature,
                                                        study)

    # calculate concordance
    concord_result_df <- calculate_RBP_KD_concordance_from_mats(test_mats_filtered,
                                                                study_mats,
                                                                signature,
                                                                study)

    # merge concordance result with linear model result
    result_df$concordance <- concord_result_df$score
    result_df$concordance_abs <- concord_result_df$score_abs
    result_df$concordance_pval <- concord_result_df$pval

    # run GSEA analysis
    # divide pos events and neg events for GSEA query
    interest_event_lib_full <- SPARKS::extract_GSEA_significant_events(test_mats)
    interest_event_lib_filtered <- SPARKS::extract_GSEA_significant_events(test_mats_filtered)

    # print(signature)
    # print(length(interest_event_lib_filtered$positive))
    # print(length(interest_event_lib_filtered$negative))

    gsea_result_full <- calculate_fgsea_score(interest_event_lib_full, study_rank)

    # skip if library length is 0
    if ((length(interest_event_lib_filtered$positive) >= 15) & (length(interest_event_lib_filtered$negative) >= 15)){
      gsea_result_filtered <- calculate_fgsea_score(interest_event_lib_filtered, study_rank)
      print(signature)
      print(gsea_result_filtered)
      # calculate rank
      gsea_pos_score <- gsea_result_filtered$pos_score
      gsea_neg_score <- gsea_result_filtered$neg_score
      gsea_pos_pval <- gsea_result_filtered$pos_pval
      gsea_neg_pval <- gsea_result_filtered$neg_pval

      gsea_score <- gsea_result_filtered$score
      gsea_score_abs <- abs(gsea_result_filtered$score)
      gsea_combined_pval <- gsea_result_filtered$pval
    } else {
      # calculate rank
      gsea_pos_score <- 0
      gsea_neg_score <- 0
      gsea_pos_pval <- 1
      gsea_neg_pval <- 1
      gsea_score <- 0
      gsea_score_abs <- 0
      gsea_combined_pval <- 1
    }

    # calculate rank
    result_df$gsea_pos_score <- gsea_pos_score
    result_df$gsea_neg_score <- gsea_neg_score
    result_df$gsea_pos_pval <- gsea_pos_pval
    result_df$gsea_neg_pval <- gsea_neg_pval
    result_df$gsea_score <- gsea_score
    result_df$gsea_score_abs <- gsea_score_abs
    result_df$gsea_combined_pval <- gsea_combined_pval

    # add tangential full information
    # calculate rank
    result_df$gsea_pos_score_full <- gsea_result_full$pos_score
    result_df$gsea_neg_score_full <- gsea_result_full$neg_score
    result_df$gsea_score_full <- gsea_result_full$score
    result_df$gsea_score_abs_full <- abs(gsea_result_full$score)
    result_df$gsea_combined_pval_full <- gsea_result_full$pval

    return(result_df)
  }, mc.cores = num_cores)

  # if running single-entry item, it gets weird, so handling it in two diff case
  # - this is a known bug in pbmcapply - issue #50
  # - we are doing "catch" approach for possible user env variations
  if (length(test_result_df_list$value) == 0) {  # this is the expected behavior
    test_result_df <- do.call(rbind, test_result_df_list)
  } else if (length(test_result_df_list$value) > 1) {  # this sometims happens because of warnings
    test_result_df <- do.call(rbind, test_result_df_list$value)
  } else {  # if single? need to be tested
    test_result_df <- test_result_df_list$value[[1]]
  }

  ### calculate RANK for downstream analysis
  if (method == "GSEA"){
    test_result_df$rank <- rank(-test_result_df$gsea_score_abs)
    test_result_df$plot_score <- test_result_df$gsea_score_abs
  } else if (method == "Pearson"){
    test_result_df$rank <- rank(-test_result_df$score_abs)
    test_result_df$plot_score <- test_result_df$score_abs
  } else if (method == "Kendall"){
    test_result_df$rank <- rank(-test_result_df$concordance_abs)
    test_result_df$plot_score <- test_result_df$concordance_abs
  }

  # add the frequency cutoff to reduce confusion downstream
  test_result_df$frequency_cutoff <- overlap_ratio_threshold
  return(test_result_df)
}



#' Add SPAKRS result for single AS signature
#'
#' @param input_mats MATS df for the study of interest
#' @param input_result SPARKS analysis result against AS signature library
#' @param test_mats MATS df for the signature
#' @param test_study signature name
#' @param kd_library If given, this entire library will be used to compute the AS events fo filter out
#' @param event_count_df If given, this list will be used to compute the AS events fo filter out. This is more memory-efficient
#' @param library_size This is 675 for now. Don't change
#'
#' @return SPARKS analysis result with new result added
#' @export
add_custom_library_to_SPARKS_test_result <- function(input_mats,
                                                     input_result,
                                                     test_mats,
                                                     test_study,
                                                     kd_library = list(),
                                                     event_count_df = data.frame(),
                                                     library_size = 675){

  # remove spl type if there is one just in case
  if ("spl_type" %in% colnames(input_result)){
    input_result$spl_type <- NULL
  }

  # remove GSEA rank if it is already there - as this will need to be recalculated
  if ("gsea_rank" %in% colnames(input_result)){
    input_result$gsea_rank <- NULL
  }

  if ((dim(event_count_df)[1] == 0) & (length(kd_library) > 0)){  # if event count is not given, update the whole library
    # add the entry to library
    added_library <- kd_library
    added_library[[test_study]] <- test_mats
  } else if ((dim(event_count_df)[1] >= 0) & (length(kd_library) == 0)) {  # if event count is given, skip the library step
    added_library <- list()
    added_library[[test_study]] <- test_mats
  } else {  # if neither or both is given - it is unclear so stop
    stop("Library input/count is missing or both present. Please check your parameter")
  }

  # run the new analysis only on the new entry
  new_test_result <- perform_SPARKS_analysis_with_overlap_filter(input_mats,
                                                                 added_library,
                                                                 study = unique(input_result$S2)[1],  # this should be one entry
                                                                 library_list = c(test_study),
                                                                 num_cores = 1,
                                                                 event_count_df = event_count_df,
                                                                 library_size = 675)
  # combine the new result
  combined_test_result <- rbind(input_result, new_test_result)

  # calculate the rank again
  combined_test_result$gsea_rank <- rank(-combined_test_result$gsea_score)

  return(combined_test_result)
}



#' Wrapper for running SPARKS for all splice types
#'
#' @param input_sparks Input SPARKS object
#' @param kd_library_all Signature library for all splice types
#' @param test_study Study name
#'
#' @return SPARKS analysis results for all splice types bound in a list
#' @export
perform_SPARKS_analysis_for_all_splice_types <- function(input_sparks,
                                                         kd_library_all,
                                                         test_study = "",
                                                         subset_group_1 = c(),
                                                         subset_group_2 = c()){

  spl_types <- c("SE", "A3SS", "A5SS")
  test_result_list <- lapply(spl_types, function(spl_type) {
    print(sprintf("Running SPARKS analysis on %s", spl_type))

    # filter raw data
    study_mats <- import_SPARKS_MATS_for_analysis(input_sparks, spl_type)

    test_result_df <- perform_SPARKS_analysis_with_overlap_filter(study_mats, kd_library_all[[spl_type]], test_study)

    return(test_result_df)
  })
  names(test_result_list) <- spl_types
  return(test_result_list)
}



##### SUBSET FUNCTIONS #####

#' Calculate new MATS df for subset of samples
#'
#' The sample column ID can be identified by PSI df stored in SPARKS object
#'
#' @param study_mats Filtered MATS df
#' @param comp_sample_id_1 Column ID for group one.
#' @param comp_sample_id_2 Column ID for group two
#' @param count_threshold Average read count for AS events. Events with read counts below this number in average is removed
#' @param min_threshold Minimum read count for AS events. Events with read counts below this number in any sample is removed
#'
#' @return subset_mats_df Filtered MATS df for the subset of the samples
#' @export
calculate_new_psi_for_subset_of_samples <- function(study_mats,
                                                    comp_sample_id_1,
                                                    comp_sample_id_2,
                                                    count_threshold = 20,
                                                    min_threshold = 20){

  # make copy to update the values
  study_mats_temp <- study_mats
  rownames(study_mats_temp) <- NULL

  comp_psi <- do.call(rbind, lapply(study_mats_temp$psi_values, function(x) as.numeric(strsplit(x, ",")[[1]][c(comp_sample_id_1, comp_sample_id_2)])))
  comp_count <- do.call(rbind, lapply(study_mats_temp$count_values, function(x) as.numeric(strsplit(x, ",")[[1]][c(comp_sample_id_1, comp_sample_id_2)])))

  # calculate PSI and count from PSI values
  comp_factor <- (study_mats_temp$skip_len / study_mats_temp$inc_len) * ((1 - comp_psi) / comp_psi)

  comp_skip <- (comp_factor / (1 + comp_factor)) * comp_count
  comp_inc <- (1 / (1 + comp_factor)) * comp_count

  # calculate new pooled count and PSI
  comp_inc_group_1 <- rowSums(comp_inc[, seq(1, length(comp_sample_id_1)), drop = FALSE], na.rm = T)
  comp_inc_group_2 <- rowSums(comp_inc[, seq(length(comp_sample_id_1) + 1, length(comp_sample_id_1) + length(comp_sample_id_2)), drop = FALSE], na.rm = T)
  comp_skip_group_1 <- rowSums(comp_skip[, seq(1, length(comp_sample_id_1)), drop = FALSE], na.rm = T)
  comp_skip_group_2 <- rowSums(comp_skip[, seq(length(comp_sample_id_1) + 1, length(comp_sample_id_1) + length(comp_sample_id_2)), drop = FALSE], na.rm = T)

  pooled_psi_group_1 <- (comp_inc_group_1 / 2) / ((comp_inc_group_1 / 2) + comp_skip_group_1)
  pooled_psi_group_2 <- (comp_inc_group_2 / 2) / ((comp_inc_group_2 / 2) + comp_skip_group_2)

  # update the values to new mats dir
  study_mats_temp$psi_values <- unlist(apply(comp_psi, 1, function(x) paste0(x, collapse=",")))
  study_mats_temp$count_values <- unlist(apply(comp_count, 1, function(x) paste0(x, collapse=",")))

  # calculate average counts
  study_mats_temp$avg_count <- rowMeans(comp_count)

  # calculate pooled PSI
  study_mats_temp$pooled_psi_1 <- pooled_psi_group_1
  study_mats_temp$pooled_psi_2 <- pooled_psi_group_2

  # calculate delta PSI and update the beta columns
  study_mats_temp$pooled_delta_psi <- pooled_psi_group_1 - pooled_psi_group_2
  study_mats_temp$beta <- pooled_psi_group_1 - pooled_psi_group_2

  # drop NA rows
  study_mats_compelte <- study_mats_temp[rowSums(is.na(study_mats_temp)) == 0, ]

  # compute pooled pval
  study_mats_temp$pooled_pval <- unlist(lapply(seq(dim(study_mats_temp)[1]), function(idx) {
    inc1 <- (comp_inc_group_1)[idx]
    inc2 <- (comp_inc_group_2)[idx]
    total1 <- (comp_inc_group_1 + comp_skip_group_1)[idx]
    total2 <- (comp_inc_group_2 + comp_skip_group_2)[idx]

    if (sum(c(inc1, inc2, total1, total2) == 0) == 0){
      a = prop.test(x = c(inc1, inc2),
                    n = c(total1, total2))
      return(a$p.value)
    } else {
      return(1)
    }
  }))

  # filter by count
  study_mats_clean <- subset(study_mats_temp, avg_count >= count_threshold )

  study_mats_clean$pval <- study_mats_clean$pooled_pval
  study_mats_clean$fdr <- p.adjust(study_mats_clean$pval, method = "BH")

  # calculate minimum_count
  study_mats_clean$min_count <- unlist(lapply(study_mats_clean$count_values,
                                              function(input_counts) min(do.call(as.numeric, strsplit(input_counts, ",")))))

  study_mats_temp <- subset(study_mats_clean, avg_count >= count_threshold & min_count >= min_threshold) %>% arrange(-beta)
  return(study_mats_temp)
}




#' Subset SPARKS object for subsets
#'
#' @param input_sparks Input SPARKS object
#' @param kd_library_all Signature library for all splice type
#' @param subset_study Study name for this subset comparison
#' @param sample_one_id Column ID for group one.
#' @param sample_two_id Column ID for group two.
#' @param num_cores Number of cores for SPARKS analysis
#' @param overlap_ratio_threshold Threshold to filter out AS events based on their occurrences in the library - default = 30%
#'
#' @return subset_SPARKS_object
#' @export
generate_subset_SPARKS_rerun <- function(input_sparks,
                                         kd_library_all,
                                         subset_study,
                                         sample_one_id,
                                         sample_two_id,
                                         num_cores = 3,
                                         overlap_ratio_threshold = 0.30){
  # make copy
  subset_sparks <- input_sparks

  # identify the index for the sample for the subset
  sample_idx_one <- sample_one_id
  sample_idx_two <- sample_two_id

  spl_types <- c("SE", "A3SS", "A5SS")
  # subset the PSI values
  dummy <- lapply(spl_types, function(spl_type){
    print(sprintf("Subsetting PSI df for %s", spl_type))
    full_psi_df <- input_sparks@psi_df[[spl_type]]
    subset_psi_df <- full_psi_df[, c(sample_idx_one, sample_idx_two)]
    # append the result for downstream analysis
    subset_sparks@psi_df[[spl_type]] <<- subset_psi_df
    return()
  })

  # subset the MATS
  dummy <- lapply(spl_types, function(spl_type){
    print(sprintf("Subsetting MATS df for %s", spl_type))
    # TODO - this has been checkd only for SE - generalize
    full_mats_df <- import_SPARKS_MATS_for_analysis(input_sparks, spl_type)
    subset_mats_df <- calculate_new_psi_for_subset_of_samples(full_mats_df,
                                                              sample_idx_one,
                                                              sample_idx_two)
    subset_mats_filtered <- subset_mats_df[rowSums(is.na(subset_mats_df)) == 0, ] %>% arrange(-beta)
    print(head(subset_mats_filtered))
    # subset_mats_filtered <- subset_mats_filtered[rev(order(subset_mats_filtered$beta)), ]
    # append the result for downstream analysis
    subset_sparks@MATS_list[[spl_type]] <<- subset_mats_filtered
  })

  print("Performing SPARKS Analysis for the subset")

  # update SPARKS analysis result - TODO - generalize?
  # subset_sparks_result <- perform_SPARKS_analysis_for_all_splice_types(subset_sparks,
  #                                                                      kd_library_all,
  #                                                                      subset_study)
  # subset_sparks@SPARKS_analysis_result <- subset_sparks_result

  subset_sparks_result_SE <- perform_SPARKS_analysis_with_overlap_filter(subset_sparks@MATS_list$SE,
                                                                         kd_library_all,
                                                                         study = subset_study,
                                                                         num_cores = num_cores,
                                                                         overlap_ratio_threshold = overlap_ratio_threshold)
  subset_sparks@SPARKS_analysis_result$SE <- subset_sparks_result_SE

  # print("Subsetting Expression df for the subset")
  # # subset the expression values
  # # calculate new index, as this may be different since it was randomized in python process
  # full_exp_df <- input_sparks@exp_df
  # exp_sample_idx_one <- which(startsWith(colnames(full_exp_df), sample_one_id))
  # exp_sample_idx_two <- which(startsWith(colnames(full_exp_df), sample_two_id))

  # subset_exp_df <- full_exp_df[, c(exp_sample_idx_one, exp_sample_idx_two)]
  # # append the result for downstream analysis
  # subset_sparks@exp_df <- subset_exp_df

  # update study
  subset_sparks@study <- subset_study

  return(subset_sparks)
}



