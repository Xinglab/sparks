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
import_SPARKS_MATS_for_analysis <- function(input_start, spl_type, count_threshold = 20){
  # process the data for analysis
  known_events <- rownames(subset(input_start@exon_annotation[[spl_type]], annotation == "Known_JC"))
  study_mats_raw <- input_start@MATS_list[[spl_type]]

  study_mats_known <- study_mats_raw[study_mats_raw$event %in% known_events, ]
  study_mats_temp <- subset(study_mats_known, avg_count >= count_threshold)
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

    test_result_df <- perform_SPARKS_analysis(study_mats, kd_library_all[[spl_type]], test_study, score_method)

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
  return(study_mats_clean)
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

