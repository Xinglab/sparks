##### FUNCTION #####
## Perform SPARKS test
#' @export
perform_SPARKS_analysis <- function(study_mats, kd_library, study, method = "GSEA", num_cores = 3){

  # divide pos events and neg events for GSEA query
  test_event_list <- extract_GSEA_significant_events(study_mats)
  test_pos_events <- test_event_list$positive
  test_neg_events <- test_event_list$negative

  # calculate correlation for signature
  test_result_df <- do.call(rbind, pbmcapply::pbmclapply(names(kd_library), function(signature){
    # test_result_df <- do.call(rbind, lapply(names(kd_library), function(signature){

    # print(signature)
    test_mats <- kd_library[[signature]]

    # calculate Rho-based enrichment score
    result_df <- calculate_RBP_KD_correlation_from_mats(test_mats, study_mats, signature, study, num_permutation = 100, beta_threshold = 0.1)

    # calculate Tau-based enrichment score
    concord_result_df <- calculate_RBP_KD_concordance_from_mats(test_mats, study_mats, signature, study)

    # merge Tau result with Rho result
    result_df$concordance <- concord_result_df$score
    result_df$concordance_abs <- concord_result_df$score_abs

    # calculate GSEA-based enrichment score
    gsea_score <- calculate_GSEA_score(test_mats, test_pos_events, test_neg_events)

    # merge GSEA result to total result
    result_df$gsea_pos_score <- gsea_score$pos_score
    result_df$gsea_neg_score <- gsea_score$neg_score
    result_df$gsea_score <- gsea_score$pos_score - gsea_score$neg_score
    result_df$gsea_score_abs <- abs(gsea_score$pos_score - gsea_score$neg_score)

    return(result_df)
  }, mc.cores = num_cores))

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

  return(test_result_df)
}


#' @export
perform_SPARKS_analysis_with_significance <- function(study_mats, kd_library, study,  method = "GSEA", num_cores = 3, num_MC = 100){

  # divide pos events and neg events for GSEA query
  n_test <- 1000

  # query null events
  event_list <- list()

  # extract all sig events in the library
  dummy <- lapply(names(kd_library), function(experiment){
    sig_events <- unlist(extract_GSEA_significant_events(kd_library[[experiment]]))
    event_list[[length(event_list) + 1]] <<- sig_events
    return()
  })

  # geenrate list and weight for them
  null_event <- unique(unlist(event_list))
  null_weight <- table(unlist(event_list))

  # calculate correlation for signature
  test_result_df <- do.call(rbind, pbmcapply::pbmclapply(names(kd_library), function(signature){
    print(signature)
    # test_result_df <- do.call(rbind, lapply(names(kd_library), function(signature){

    # print(signature)
    test_mats <- kd_library[[signature]]

    # calculate Rho-based enrichment score
    result_df <- calculate_RBP_KD_correlation_from_mats(test_mats, study_mats, signature, study, num_permutation = 100, beta_threshold = 0.1)

    # calculate Tau-based enrichment score
    concord_result_df <- calculate_RBP_KD_concordance_from_mats(test_mats, study_mats, signature, study)

    # merge Tau result with Rho result
    result_df$concordance <- concord_result_df$score
    result_df$concordance_abs <- concord_result_df$score_abs

    # calculate GSEA-based enrichment score
    # flip the direciton
    test_mats_event_list <- extract_GSEA_significant_events(test_mats)
    test_mats_pos <- test_mats_event_list$positive
    test_mats_neg <- test_mats_event_list$negative

    # run the analysis for N times
    score_dist <- do.call(rbind, lapply(seq(1, n_test), function(ind){
      gsea_score_pos_neg <- calculate_GSEA_score_sampling(study_mats, test_mats_pos, test_mats_neg, sampling_size = num_MC)
      gsea_score_non_sig <- calculate_GSEA_score_sampling_1d_weight(study_mats, null_event, null_weight, sampling_size = num_MC)

      score_df <- data.frame(positive = gsea_score_pos_neg$pos_score,
                             negative = gsea_score_pos_neg$neg_score,
                             null = gsea_score_non_sig)
    }))

    pos_sig_dist <- score_dist$positive
    neg_sig_dist <- score_dist$negative
    nonsig_dist <- score_dist$null

    ## calculate permutation p-value
    # calculate null values
    null_upper_score <- quantile(nonsig_dist, c(0.95))
    null_lower_score <- quantile(nonsig_dist, c(0.05))
    null_mean_score <- mean(nonsig_dist)

    # calculate both directions for pos and neg
    pos_upper_count <- sum(pos_sig_dist > null_upper_score)
    pos_lower_count <- sum(pos_sig_dist < null_lower_score)
    pos_mean_score <- mean(pos_sig_dist)

    neg_upper_count <- sum(neg_sig_dist > null_upper_score)
    neg_lower_count <- sum(neg_sig_dist < null_lower_score)
    neg_mean_score <- mean(neg_sig_dist)

    if(pos_upper_count + neg_lower_count > pos_lower_count + neg_upper_count){ # if positive enrichment
      p_perm_pos <- max(1, (n_test - pos_upper_count)) / n_test # technically p-value lower bound is 1/N
      p_perm_neg <- max(1, (n_test - neg_lower_count)) / n_test
    } else { # if negative enrichment
      p_perm_pos <- max(1, (n_test - pos_lower_count)) / n_test
      p_perm_neg <- max(1, (n_test - neg_upper_count)) / n_test
    }
    p_perm <- metap::sumlog(c(p_perm_pos, p_perm_neg))$p

    # calculate total score
    total_score <- pos_mean_score - neg_mean_score



    # merge GSEA result to total result
    result_df$gsea_pos_score <- pos_mean_score
    result_df$gsea_neg_score <- neg_mean_score
    result_df$gsea_score <- total_score
    result_df$gsea_score_abs <- abs(total_score)
    result_df$gsea_null_upper <- null_upper_score
    result_df$gsea_null_lower <- null_lower_score
    result_df$gsea_pos_upper <- pos_upper_count
    result_df$gsea_pos_lower <- pos_lower_count
    result_df$gsea_neg_upper <- neg_upper_count
    result_df$gsea_neg_lower <- neg_lower_count
    result_df$gsea_combined_pval <- p_perm


    return(result_df)
  }, mc.cores = num_cores))

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

  return(test_result_df)
}


calculate_GSEA_score_sampling <- function(input_signature_list, test_pos_events, test_neg_events, sampling_size = 100){
  pos_overlap_list <- test_pos_events[test_pos_events %in% input_signature_list$event]
  pos_overlap_list_sample <- sample(pos_overlap_list, min(sampling_size, length(pos_overlap_list)))
  if (length(pos_overlap_list) > 0){
    pos_result <- GSEA.EnrichmentScore(input_signature_list$event, pos_overlap_list_sample, weighted.score.type = 0)
    pos_stat <- pos_result$ES
  } else {
    pos_stat <- 0
  }

  neg_overlap_list <- test_neg_events[test_neg_events %in% input_signature_list$event]
  neg_overlap_list_sample <- sample(neg_overlap_list, min(sampling_size, length(neg_overlap_list)))

  if (length(neg_overlap_list) > 0){
    neg_result <- GSEA.EnrichmentScore(input_signature_list$event, neg_overlap_list_sample, weighted.score.type = 0)
    neg_stat <- neg_result$ES
  } else {
    neg_stat <- 0
  }

  result_df <- data.frame(pos_score = pos_stat,
                          neg_score = neg_stat)
  return(result_df)
}


calculate_GSEA_score_sampling_1d_weight <- function(input_signature_list, non_sig_events, weight, sampling_size = 100){
  overlap_list <- non_sig_events[non_sig_events %in% input_signature_list$event]

  tot_overlap_list_sample <- sample(overlap_list, min(sampling_size, length(overlap_list)), prob = weight[overlap_list])

  if (length(tot_overlap_list_sample) > 0){
    neg_result <- GSEA.EnrichmentScore(input_signature_list$event, tot_overlap_list_sample, weighted.score.type = 0)
    neg_stat <- neg_result$ES
  } else {
    neg_stat <- 0
  }

  score <- neg_stat
  return(score)
}




#' @export
import_SPARKS_MATS_for_analysis <- function(input_start, spl_type, count_threshold = 20){
  # check if the event name is already processed
  if(length(strsplit(input_start@MATS_list[[spl_type]][1, 'event'], ":")[[1]]) != 8){  # if not processed
    # process the data for analysis
    known_events <- rownames(subset(input_start@exon_annotation[[spl_type]], annotation == "Known_JC"))
    study_mats_raw <- input_start@MATS_list[[spl_type]]

    study_mats_known <- study_mats_raw[study_mats_raw$event %in% known_events, ]

    # update the event to short form
    study_mats_known$event <- unlist(lapply(study_mats_known$event,
                                            function(x) rewrite_event_coordinates(x)))
  } else {  # means the MATS is already processed, so no need to anything
    study_mats_known <- input_start@MATS_list[[spl_type]]
  }

  # calculate minimum_count
  study_mats_known$min_count <- unlist(lapply(study_mats_known$count_values,
                                              function(input_counts) min(do.call(as.numeric, strsplit(input_counts, ",")))))

  # filter by minimum count
  study_mats_temp <- subset(study_mats_known,
                            avg_count >= count_threshold &
                              min_count >= count_threshold) %>% arrange(-beta)
  return(study_mats_temp)
}


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


#' @export
perform_SPARKS_analysis_for_all_splice_types <- function(input_start, kd_library_all,
                                                         test_study = "", score_method = "GSEA",
                                                         subset_group_1 = c(), subset_group_2 = c()){

  spl_types <- c("SE", "A3SS", "A5SS")

  test_result_list <- lapply(spl_types, function(spl_type) {
    # filter raw data
    study_mats <- import_SPARKS_MATS_for_analysis(input_start, spl_type)

    if (length(subset_group_1) > 0 & length(subset_group_2) > 0){
      study_mats <- calculate_new_psi_for_subset_of_samples(study_mats, subset_group_1, subset_group_2)
    }

    test_result_df <- perform_SPARKS_analysis_with_overlap_filter(study_mats, kd_library_all[[spl_type]], test_study, score_method)

    # store data
    test_result_df$spl_type <- spl_type

    return(test_result_df)
  })
  names(test_result_list) <- spl_types
  return(test_result_list)
}


#' @export
generate_SPARKS_top10_plots <- function(test_result_list){

  spl_types <- c("SE", "A3SS", "A5SS", "RI")
  ### generate plot data
  # subset top 10

  test_plot_list <- lapply(spl_types, function(spl_type) {
    test_result_df <- test_result_list[[spl_type]]
    test_top10_result_df <- subset(test_result_df, rank <= 10)

    test_top10_result_df$S1 <- factor(test_top10_result_df$S1,
                                      level = test_top10_result_df$S1[order(test_top10_result_df$rank)])
    p <- ggplot(test_top10_result_df,
                aes(x = rank,
                    y = plot_score)) +
      geom_point() +
      ggrepel::geom_label_repel(data = subset(test_top10_result_df, rank < 6),
                                hjust = -0.1, vjust = -0.1, aes(label = S1)) +
      labs(x = 'Rank based on Correlation Score from ENCODE RBP KD \n(Sorted by score)',
           y = "Correlation Score",
           subtitle = sprintf("Top 10 RBP KD experiments based on Correlation Score\nSplicing Type: %s", spl_type)) +
      scale_x_continuous(breaks = seq(10))
    return(p)
  })
  return(test_plot_list)
}


#' @export
generate_SPARKS_sorted_barplot_with_highlight <- function(test_result_list, spl_type, rbp_of_interest){

  # extract relevant score
  # TODO - score scheme should be adaptable to different method
  test_result_df <- test_result_list[[spl_type]]
  test_result_df$plot_scores <- test_result_df$gsea_score

  # sort by plot score
  test_result_df$S1 <- factor(test_result_df$S1, levels = test_result_df[order(test_result_df$plot_scores), ]$S1)

  # annotate RBP and flag for highlight
  test_result_df$RBP <- unlist(lapply(as.character(test_result_df$S1), function(x) strsplit(x, "_")[[1]][2]))
  test_result_df$interest <- ifelse(test_result_df$RBP %in% c(rbp_of_interest), "interest", NA)

  # generate plot
  p <- ggplot(test_result_df, aes(x = S1, y = plot_scores, fill = interest)) +
    geom_bar(stat = 'identity', width = 1) +
    ggrepel::geom_label_repel(data = subset(test_result_df, interest == "interest"),
                              aes(label = S1),
                              vjust = 1, hjust = -0.1,
                              min.segment.length = 0) +
    scale_fill_manual(values = list('interest' = "indianred3"),
                      na.value = 'grey75') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    labs(y = "Correlation Score",
         x = "ENCODE RBP KD/KO experiments\n(Sorted by Correlation Score)")
  return(p)
}


#' @export
## Subset data if necessary
calculate_new_psi_for_subset_of_samples <- function(study_mats, comp_sample_id_1, comp_sample_id_2, count_threshold = 20){

  # make copy to update the values
  study_mats_temp <- study_mats
  rownames(study_mats_temp) <- NULL

  comp_psi <- do.call(rbind, lapply(study_mats_temp$psi_values, function(x) as.numeric(strsplit(x, ",")[[1]][c(comp_sample_id_1, comp_sample_id_2)])))
  comp_count <- do.call(rbind, lapply(study_mats_temp$count_values, function(x) as.numeric(strsplit(x, ",")[[1]][c(comp_sample_id_1, comp_sample_id_2)])))

  # calculate PSI and count from PSI values
  comp_factor <- (study_mats_temp$skip_len / study_mats_temp$inc_len) * ((1 - comp_psi) / comp_psi)

  comp_skip <- (comp_factor / (1 + comp_factor)) * comp_count
  comp_inc <- (1 / (1 + comp_factor)) * comp_count


  # comp_inc <- comp_psi * comp_count * 2 / (1 + comp_psi)
  # comp_skip <- comp_count - comp_inc

  # # calculate loss
  # comp_psi_new = (comp_inc / 2) / ((comp_inc / 2) + comp_skip)
  # max(comp_psi - comp_psi_new, na.rm = T)

  # calculate new pulled count and PSI
  comp_inc_group_1 <- rowSums(comp_inc[, seq(1, length(comp_sample_id_1)), drop = FALSE], na.rm = T)
  comp_inc_group_2 <- rowSums(comp_inc[, seq(length(comp_sample_id_1) + 1, length(comp_sample_id_1) + length(comp_sample_id_2)), drop = FALSE], na.rm = T)
  comp_skip_group_1 <- rowSums(comp_skip[, seq(1, length(comp_sample_id_1)), drop = FALSE], na.rm = T)
  comp_skip_group_2 <- rowSums(comp_skip[, seq(length(comp_sample_id_1) + 1, length(comp_sample_id_1) + length(comp_sample_id_2)), drop = FALSE], na.rm = T)

  pulled_psi_group_1 <- (comp_inc_group_1 / 2) / ((comp_inc_group_1 / 2) + comp_skip_group_1)
  pulled_psi_group_2 <- (comp_inc_group_2 / 2) / ((comp_inc_group_2 / 2) + comp_skip_group_2)

  # update the values to new mats dir
  study_mats_temp$psi_values <- unlist(apply(comp_psi, 1, function(x) paste0(x, collapse=",")))
  study_mats_temp$count_values <- unlist(apply(comp_count, 1, function(x) paste0(x, collapse=",")))

  study_mats_temp$avg_count <- rowMeans(comp_count)

  study_mats_temp$pulled_psi_1 <- pulled_psi_group_1
  study_mats_temp$pulled_psi_2 <- pulled_psi_group_2

  study_mats_temp$pulled_delta_psi <- pulled_psi_group_1 - pulled_psi_group_2
  study_mats_temp$beta <- pulled_psi_group_1 - pulled_psi_group_2

  # drop NA rows
  study_mats_compelte <- study_mats_temp[rowSums(is.na(study_mats_temp)) == 0, ]

  study_mats_temp$pulled_pval <- unlist(lapply(seq(dim(study_mats_temp)[1]), function(idx) {
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

  study_mats_clean$pval <- study_mats_clean$pulled_pval
  study_mats_clean$fdr <- p.adjust(study_mats_clean$pval, method = "BH")

  # calculate minimum_count
  study_mats_clean$min_count <- unlist(lapply(study_mats_clean$count_values,
                                              function(input_counts) min(do.call(as.numeric, strsplit(input_counts, ",")))))

  study_mats_temp <- subset(study_mats_clean, avg_count >= count_threshold & min_count >= count_threshold) %>% arrange(-beta)
  return(study_mats_temp)
}




#' @export
generate_SPARKS_strip_plot_with_highlight <- function(test_result_list, spl_type, rbp_of_interest){

  # extract relevant score
  # TODO - score scheme should be adaptable to different method
  test_result_df <- test_result_list[[spl_type]]

  # sort by plot score
  test_result_df$plot_scores <- test_result_df$gsea_score
  test_result_df$S1 <- factor(test_result_df$S1, levels = test_result_df[order(test_result_df$plot_scores), ]$S1)
  test_result_df$plot_rank <- rank(-test_result_df$plot_scores)

  # annotate RBP and flag for highlight
  test_result_df$RBP <- unlist(lapply(as.character(test_result_df$S1), function(x) strsplit(x, "_")[[1]][2]))
  test_result_df$interest <- ifelse(test_result_df$RBP %in% c(rbp_of_interest), "interest", NA)
  test_result_df[!is.na(test_result_df$interest), ]$plot_scores <- NA  # Using NA to for manual setting

  # set names for ranking
  test_result_df$annotation <- unlist(apply(test_result_df, 1, function(x) {
    S1 <- x["S1"]
    plot_rank <- round(as.numeric(x["plot_rank"]))
    sprintf("%s - %s %s %s",
            plot_rank,
            strsplit(as.character(S1), "_")[[1]][1],
            strsplit(as.character(S1), "_")[[1]][2],
            strsplit(as.character(S1), "_")[[1]][3])
  }))



  # generate plot
  p <- ggplot(test_result_df, aes(y = S1)) +
    geom_tile(aes(x = 0.05,
                  fill = plot_scores),
              width = 0.1) +
    scale_fill_distiller(palette = "PiYG", limits = c(-max(abs(test_result_df$gsea_score)),
                                                      max(abs(test_result_df$gsea_score))),
                         na.value = "black")  +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 0.3)) +
    scale_y_discrete(expand = c(0, 0)) +
    ggrepel::geom_label_repel(data = subset(test_result_df, interest == "interest"),
                              aes(x = 0.1 - 0.000005, # this is to remove slight disconnect between the bar and label line
                                  label = annotation,
                                  fill = gsea_score),
                              # vjust = 1,
                              # hjust = 0.5,
                              min.segment.length = 0,
                              nudge_y = -20,
                              nudge_x = 0.05,
                              force = 10,
                              xlim = c(0.1, 0.3),
                              ylim = c(5, dim(test_result_df)[1] - 5),
                              label.r = 0,
                              label.padding = 0.3,
                              family = "Arial",
                              label.size = NA,
                              size = 3) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.direction = "vertical",
          axis.title = element_blank()) +
    labs(x = "Enrichment Score\nFrom the Library",
         y = spl_type)

  return(p)
}

#' @export
generate_splicing_summary_plot <- function(study_mats_list){

  test_sig_df  <- do.call(rbind, lapply(names(study_mats_list), function(spl_type){
    test_mats <- study_mats_list[[spl_type]]
    test_sig_subset <- subset(test_mats, abs(beta) > 0.1)
    test_sig_subset$spl_type <- spl_type

    return(test_sig_subset)
  }))


  test_sig_count_df  <- data.frame(table(data.frame(ifelse(test_sig_df $beta > 0, "pos", "neg"), test_sig_df $spl_type)))
  colnames(test_sig_count_df ) <- c("direction", "spl_type", "freq")
  test_sig_count_df $loc <- ifelse(test_sig_count_df $direction == "neg", max(test_sig_df $beta) + 0.1, min(test_sig_df $beta) - 0.1)

  p <- ggplot(test_sig_df) +
    geom_violin(aes(x = spl_type, fill = ifelse(beta > 0, "pos", "neg"),
                    y = beta * -1), scale = "width" ) +
    # geom_jitter(alpha = 0.1)
    geom_point(aes(x = spl_type, fill = ifelse(beta > 0, "pos", "neg"),
                   y = beta* -1),
               position=position_jitterdodge(), alpha = 0.1) +
    geom_label(data = test_sig_count_df , aes(label = freq, x = spl_type, y = loc, color = direction), fill = "white") +
    labs(x = "Splicing Type",
         y = "Percent Spliced In") +
    theme(legend.position = "none")
  return(p)
}


#' @export
prepare_library_data_for_splicing_summary_plot <- function(kd_library_all, study){

  output_list <- list()
  spl_types <- c("SE", "A3SS", "A5SS", "RI")

  dummy <- lapply(spl_types, function(spl_type){
    study_mats <- kd_library_all[[spl_type]][[study]]
    output_list[[spl_type]] <<- study_mats
  })
  return(output_list)
}

#' @export
add_custom_study_to_library <- function(kd_library_all,
                                        input_file,
                                        input_study,
                                        subset_group_1 = c(),
                                        subset_group_2 = c()){
  # read SPARKS data
  input_start <- readRDS(input_file)

  spl_types <- c("SE", "A3SS", "A5SS")

  dummy <- lapply(spl_types, function(spl_type){
    print(sprintf("Importing %s for %s",spl_type, input_study))

    input_mats_filtered <- import_SPARKS_MATS_for_analysis(input_start, spl_type = spl_type, count_threshold = 20)

    # subset option if input necessitates it
    if (length(subset_group_1) > 0 | length(subset_group_2) > 0){  # if subset option is used
      # validate the length for both
      if (length(subset_group_1) == 0 | length(subset_group_2) == 0){  # if one group is not specified
        stop("One of the subset group is undefined - please check your input")
      }
      # perform subset operation
      input_mats_filtered <- calculate_new_psi_for_subset_of_samples(input_mats_filtered,
                                                                     subset_group_1,
                                                                     subset_group_2)
    }

    input_mats_sorted <- input_mats_filtered[rev(order(input_mats_filtered$beta)), ]


    concordant_mats <- subset(input_mats_sorted, abs(beta - pulled_delta_psi) < 0.1)

    # change beta
    concordant_mats$mean_beta <- concordant_mats$beta
    concordant_mats$beta <- concordant_mats$pulled_delta_psi


    kd_library_all[[spl_type]][[input_study]] <<- concordant_mats
    # discordant_mats_list[[study]] <<- discordant_mats
    return()
  })

  return(kd_library_all)
}

#' @export
import_custom_study_as_SPARKS_library <- function(input_SPARKS_file,
                                                  input_study,
                                                  subset_group_1 = c(),
                                                  subset_group_2 = c()){
  # read SPARKS data
  input_start <- readRDS(input_SPARKS_file)

  spl_types <- c("SE", "A3SS", "A5SS")

  concordant_mats_list <- list()

  # make slots
  concordant_mats_list[['SE']] <- list()
  concordant_mats_list[['A3SS']] <- list()
  concordant_mats_list[['A5SS']] <- list()


  dummy <- lapply(spl_types, function(spl_type){
    print(sprintf("Importing %s for %s", spl_type, input_study))

    input_mats_filtered <- import_SPARKS_MATS_for_analysis(input_start, spl_type = spl_type, count_threshold = 20)

    # subset option if input necessitates it
    if (length(subset_group_1) > 0 | length(subset_group_2) > 0){  # if subset option is used
      # validate the length for both
      if (length(subset_group_1) == 0 | length(subset_group_2) == 0){  # if one group is not specified
        stop("One of the subset group is undefined - please check your input")
      }
      # perform subset operation
      input_mats_filtered <- calculate_new_psi_for_subset_of_samples(input_mats_filtered,
                                                                     subset_group_1,
                                                                     subset_group_2)
    }

    input_mats_sorted <- input_mats_filtered[rev(order(input_mats_filtered$beta)), ]


    concordant_mats <- subset(input_mats_sorted, abs(beta - pulled_delta_psi) < 0.1)

    # change beta
    concordant_mats$mean_beta <- concordant_mats$beta
    concordant_mats$beta <- concordant_mats$pulled_delta_psi


    concordant_mats_list[[spl_type]][[input_study]] <<- concordant_mats
    # discordant_mats_list[[study]] <<- discordant_mats
    return()
  })

  return(concordant_mats_list)
}

#' @export
generate_SPARKS_result_with_custom_library_results <- function(input_sparks,
                                                               custom_test_result){

  input_test_result <- input_sparks@SPARKS_analysis_result
  combined_test_result <- add_custom_library_result_to_SPARKS_result(input_test_result,
                                                                     custom_test_result)
  return(combined_test_result)
}

#' @export
add_custom_library_result_to_SPARKS_result <- function(input_test_result,
                                                       custom_test_result){
  # add this to the original library data
  combined_test_result <- list()
  dummy <- lapply(spl_types, function(spl_type){
    # load result for each type
    custom_test_result_spl_type <- custom_test_result[[spl_type]]
    original_test_result_spl_type <- input_test_result[[spl_type]]

    # combine the result
    combined_result_spl_type <- rbind(custom_test_result_spl_type,
                                      original_test_result_spl_type)

    # recalculate rank
    combined_result_spl_type$rank <- ceiling(rank(-combined_result_spl_type$gsea_score))  # TODO - change this code for other methods?

    # update the dataframe
    combined_test_result[[spl_type]] <<- combined_result_spl_type
    return()
  })
  return(combined_test_result)
}

#' @export
merge_custom_SPARKS_libraries <- function(custom_library_list){

  # read SPARKS data
  spl_types <- c("SE", "A3SS", "A5SS")

  merged_mats_list <- list()

  # make slots
  merged_mats_list[['SE']] <- list()
  merged_mats_list[['A3SS']] <- list()
  merged_mats_list[['A5SS']] <- list()

  input_study_list <- unlist(lapply(custom_library_list, function(x) names(x[['SE']])))

  dummy <- lapply(spl_types, function(spl_type){
    dummy2 <- lapply(seq(length(custom_library_list)), function(idx){  # access by index for easier processing
      input_study <- input_study_list[idx]
      print(sprintf("Merging %s for %s", spl_type, input_study))

      # access the correct study mats
      library_mats <- custom_library_list[[idx]][[spl_type]][[input_study]]

      merged_mats_list[[spl_type]][[input_study]] <<- library_mats

      return()
    })
    return()
  })

  return(merged_mats_list)
}



#' @export
append_SPARKS_mats_for_all_splice_type <- function(mats_list,
                                                   input_start,
                                                   input_study){
  # make the slot for input
  mats_list[[input_study]] <- list()

  dummy <- lapply(spl_types, function(spl_type){
    input_study_mats <- import_SPARKS_MATS_for_analysis(input_start, spl_type)

    # append the result for downstream analysis
    mats_list[[input_study]][[spl_type]] <<- input_study_mats
  })

  return(mats_list)
}


#' @export
append_SPARKS_result_for_all_splice_type <- function(result_list,
                                                     input_start,
                                                     input_study){
  # make the slot for input
  result_list[[input_study]] <- list()

  input_result_list <- input_start@SPARKS_analysis_result

  # append the result for downstream analysis
  result_list[[input_study]] <- input_result_list


  return(result_list)
}



#' @export
append_SPARKS_expression <- function(expression_list,
                                     input_start,
                                     input_study){
  # make the slot for input
  expression_list[[input_study]] <- list()

  input_expression_list <- input_start@exp_df

  # append the expression for downstream analysis
  expression_list[[input_study]] <- input_expression_list


  return(expression_list)
}



#' @export
generate_subset_SPARKS <- function(input_start,
                                   kd_library_all,
                                   subset_study,
                                   sample_one_id,
                                   sample_two_id){
  # make copy
  subset_start <- input_start

  # identify the index for the sample for the subset
  sample_idx_one <- which(startsWith(rownames(input_start@sample_chart), sample_one_id))
  sample_idx_two <- which(startsWith(rownames(input_start@sample_chart), sample_two_id))

  # subset the PSI values
  dummy <- lapply(spl_types, function(spl_type){
    print(sprintf("Subsetting PSI df for %s", spl_type))
    full_psi_df <- input_start@psi_df[[spl_type]]
    subset_psi_df <- full_psi_df[, c(sample_idx_one, sample_idx_two)]
    # append the result for downstream analysis
    subset_start@psi_df[[spl_type]] <<- subset_psi_df
    return()
  })

  # subset the MATS
  dummy <- lapply(spl_types, function(spl_type){
    print(sprintf("Subsetting MATS df for %s", spl_type))

    full_mats_df <- import_SPARKS_MATS_for_analysis(input_start, spl_type)
    subset_mats_df <- calculate_new_psi_for_subset_of_samples(full_mats_df,
                                                              sample_idx_one,
                                                              sample_idx_two)
    # append the result for downstream analysis
    subset_start@MATS_list[[spl_type]] <<- subset_mats_df
  })

  print("Performing SPARKS Analysis for the subset")

  # update SPARKS analysis result
  subset_sparks_result <- perform_SPARKS_analysis_for_all_splice_types(subset_start,
                                                                       kd_library_all,
                                                                       subset_study)
  subset_start@SPARKS_analysis_result <- subset_sparks_result

  print("Subsetting Expression df for the subset")
  # subset the expression values
  # calculate new index, as this may be different since it was randomized in python process
  full_exp_df <- input_start@exp_df
  exp_sample_idx_one <- which(startsWith(colnames(full_exp_df), sample_one_id))
  exp_sample_idx_two <- which(startsWith(colnames(full_exp_df), sample_two_id))

  subset_exp_df <- full_exp_df[, c(exp_sample_idx_one, exp_sample_idx_two)]
  # append the result for downstream analysis
  subset_start@exp_df <- subset_exp_df

  # update study
  subset_start@study <- subset_study

  return(subset_start)
}


#' @export
#'
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
    #   # we only make the p-value 1
    #   combined_pval <- 1
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



#' @export
perform_SPARKS_analysis_with_overlap_filter <- function(study_mats,
                                                        kd_library,
                                                        study,
                                                        method = "GSEA",
                                                        num_cores = 3,
                                                        overlap_ratio_threshold = 0.20,
                                                        library_list = c()){
  # query null events
  event_list <- list()

  # extract all sig events in the library
  dummy <- lapply(names(kd_library), function(experiment){
    sig_events <- unlist(extract_GSEA_significant_events(kd_library[[experiment]]))
    event_list[[length(event_list) + 1]] <<- sig_events
    return()
  })

  # generate list and weight for them
  null_event <- unique(unlist(event_list))
  null_weight <- table(unlist(event_list))

  # define uninformative events
  uninformative_events <- names(null_weight)[null_weight > length(kd_library) * overlap_ratio_threshold]

  study_mats_clean <- study_mats[, c("event", "pulled_delta_psi")]
  study_rank <- study_mats_clean$pulled_delta_psi
  names(study_rank) <- study_mats_clean$event

  # limit the scope of kd library if given
  if(length(library_list) > 0){
    kd_library <- kd_library[library_list]
  }

  # calculate correlation for signature
  test_result_df <- do.call(rbind, pbmcapply::pbmclapply(names(kd_library), function(signature){

    print(signature)
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

    # print("AAAA")
    # run GSEA analysis
    # divide pos events and neg events for GSEA query
    interest_event_lib_full <- SPARKS::extract_GSEA_significant_events(test_mats)
    interest_event_lib_filtered <- SPARKS::extract_GSEA_significant_events(test_mats_filtered)

    print(signature)
    print(length(interest_event_lib_filtered$positive))
    print(length(interest_event_lib_filtered$negative))

    # print(library_sig_interest)
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


    ### PERFORM permutation
    ## Permutation is disabled because it doesn't really add much value
    ## compared to traditional BH FDR correction
    ## but it takes a lot of computational resources
    # # get size
    # size_pos <- length(interest_event_lib_filtered$positive)
    # size_neg <- length(interest_event_lib_filtered$negative)
    #
    # pval_list <- lapply(seq(n_test), function(idx){
    #
    #   # generate size matching null
    #   pos_null_events <- sample(names(study_rank), size_pos)
    #   neg_null_events <- sample(names(study_rank), size_neg)
    #
    #   library_null <- list()
    #   library_null$positive <- pos_null_events
    #   library_null$negative <- neg_null_events
    #
    #   # run gsea
    #   gsea_result_null <- calculate_fgsea_score(library_null, study_rank)
    #
    #   combined_pval <- gsea_result_null$pval
    #
    #   return(combined_pval)
    # })
    #
    # perm_pval <- min(unlist(pval_list))




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

    # add perm pval
    # result_df$perm_pval <- perm_pval
    # print(result_df)

    return(result_df)
  }, mc.cores = num_cores))

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

  return(test_result_df)
}


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
    } else if (elements[[4]] == "-"){  # if negavite
      new_event <- paste0(c(elements[2:4], elements[9], elements[6], elements[8], elements[11]), collapse = ":")
    } else {
      stop()
    }
  } else if (spl_type == "A5SS"){
    if (elements[[4]] == '+'){  # if positive
      new_event <- paste0(c(elements[2:4], elements[8], elements[6], elements[9], elements[11]), collapse = ":")
    } else if (elements[[4]] == "-"){  # if negavite
      new_event <- paste0(c(elements[2:4], elements[7], elements[5], elements[10], elements[11]), collapse = ":")
    } else {
      stop()
    }
  }
  return(new_event)
}


#' @export
import_SPARKS_MATS_for_rerun <- function(input_spark, spl_type = "SE"){
  # import study mats
  study_mats <- import_SPARKS_MATS_for_analysis(input_spark, spl_type)


  # check if the event name is already processed
  if(length(strsplit(study_mats[1, 'event'], ":")[[1]]) != 8){  # if not processed
    # clean the coordinates
    study_mats$event <- unlist(lapply(study_mats$event,
                                      function(x) rewrite_event_coordinates(x)))
  }

  # sort study mats to use as profile
  # study_mats_sorted <- study_mats[order(-study_mats$beta), ]
  study_mats_sorted <- study_mats %>% arrange(-beta)

  return(study_mats_sorted)
}


#' @export
generate_subset_SPARKS_rerun <- function(input_sparks,
                                         kd_library_all,
                                         subset_study,
                                         sample_one_id,
                                         sample_two_id,
                                         num_cores = 3,
                                         num_MC = 300,
                                         overlap_ratio_threshold = 0.20){
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
    full_mats_df <- import_SPARKS_MATS_for_rerun(input_sparks, spl_type)
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

