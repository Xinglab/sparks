#' Generate scatter plot with linear fit statistics between two AS changes
#'
#' @param s1_mats_df Filtered MATS df for study 1
#' @param s2_mats_df Filtered MATS df for study 2
#' @param study_1 Name for study 1
#' @param study_2 Name for study 2
#' @param beta_threshold Minimum delPSI threshold. AS events below this for either studies are discarded.
#' @param slope If TRUE - show slope of the linear fit
#'
#' @return Scatter plot of AS changes betwene those two studies
#' @export
generate_RBP_KD_correlation_scatter_plot <- function(s1_mats_df,
                                                     s2_mats_df,
                                                     study_1,
                                                     study_2,
                                                     beta_threshold = 0.1,
                                                     slope = FALSE){

  # merge each mats by event
  merged_mats_df <- merge(s1_mats_df, s2_mats_df, by = "event")

  delta_psi_df_prefilter <- merged_mats_df[, c("beta.x", "beta.y")]
  rownames(delta_psi_df_prefilter) <- merged_mats_df$event
  colnames(delta_psi_df_prefilter) <- c("S1", "S2")

  # calculate number of total observations
  num_complete <- sum(rowSums(!is.na(delta_psi_df_prefilter)) == 2 )
  num_S1_total <- sum(!is.na(delta_psi_df_prefilter$S1))
  num_S2_total <- sum(!is.na(delta_psi_df_prefilter$S2))

  # calculate significant events regardless of common events
  num_S1_sig <- sum(abs(delta_psi_df_prefilter$S1) > beta_threshold)
  num_S2_sig <- sum(abs(delta_psi_df_prefilter$S2) > beta_threshold)

  # keep only complete observations
  delta_psi_df <- delta_psi_df_prefilter[rowSums(is.na(delta_psi_df_prefilter)) == 0, ]

  # subset only significant ones
  # only keep events with delta psi above threshold in either studies for count data
  psi_df_any_sig <- subset(delta_psi_df, abs(S1) > beta_threshold | abs(S2) > beta_threshold)
  # use both sig events for linear model
  psi_df_both_sig <- subset(delta_psi_df, abs(S1) > beta_threshold & abs(S2) > beta_threshold)

  # calculate summary data to find beta direction
  num_pos_pos <- dim(subset(psi_df_any_sig, S1 > beta_threshold & S2 > beta_threshold))[1]
  num_pos_neg <- dim(subset(psi_df_any_sig, S1 > beta_threshold & S2 < -beta_threshold))[1]
  num_neg_pos <- dim(subset(psi_df_any_sig, S1 < -beta_threshold & S2 > beta_threshold))[1]
  num_neg_neg <- dim(subset(psi_df_any_sig, S1 < -beta_threshold & S2 < -beta_threshold))[1]

  # calculate neutral ones for filtering as well
  num_not_pos <- dim(subset(psi_df_any_sig, S1 < beta_threshold & S1 > -beta_threshold & S2 > beta_threshold))[1]
  num_not_neg <- dim(subset(psi_df_any_sig, S1 < beta_threshold & S1 > -beta_threshold & S2 < -beta_threshold))[1]
  num_pos_not <- dim(subset(psi_df_any_sig, S1 > beta_threshold & S2 < beta_threshold & S2 > -beta_threshold))[1]
  num_neg_not <- dim(subset(psi_df_any_sig, S1 < -beta_threshold & S2 < beta_threshold & S2 > -beta_threshold))[1]


  # run linear model and get relevant statistics
  tryCatch({
    # lm_result <- lm(psi_df_any_sig$S2 ~ psi_df_any_sig$S1 + 0)
    lm_result <- lm(psi_df_both_sig$S2 ~ psi_df_both_sig$S1 + 0)

    rsq <- summary(lm_result)$r.squared
    correlation_value <- as.numeric(lm_result$coefficients)
    summary_stats <- summary(lm_result)$fstatistic
    pval <- pf(summary_stats[1], summary_stats[2], summary_stats[3], lower.tail=F)
  },
  error = function(e){
    rsq <<- 0
    correlation_value <<- 0
    pval <<- 1
  })

  # # diagnostic plot for dev
  # ggplot(delta_psi_df, aes(x = S1, y = S2)) + geom_point()
  #

  # calculate effective number
  num_eff <- (num_pos_pos + num_neg_neg) - (num_neg_pos + num_pos_neg)
  num_total <- (num_pos_pos + num_neg_neg) + (num_neg_pos + num_pos_neg) + (num_not_pos + num_not_neg) + (num_pos_not + num_neg_not)

  # calculate score
  score <- rsq * num_eff / num_total
  score_abs <- abs(score)

  ## generate plot
  # define coordinates
  # manually set 4 ticks in each direction
  max_coord <- ceiling(max(min(psi_df_any_sig$S1, psi_df_any_sig$S2), max(psi_df_any_sig$S1, psi_df_any_sig$S2)) * 4) / 4
  # min_coord = -max(min(psi_df_any_sig$S1, psi_df_any_sig$S2), max(psi_df_any_sig$S1, psi_df_any_sig$S2))
  min_coord <- -1 * max_coord


  p <- ggplot(psi_df_any_sig, aes(x = S1, y = S2)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = beta_threshold, color = "indianred2", linetype = "dashed") +
    geom_hline(yintercept = -beta_threshold, color = "indianred2", linetype = "dashed") +
    geom_vline(xintercept = beta_threshold, color = "indianred2", linetype = "dashed") +
    geom_vline(xintercept = -beta_threshold, color = "indianred2", linetype = "dashed") +
    scale_x_continuous(limits = c(min_coord, max_coord), breaks = seq(-ceiling(max_coord*4), ceiling(max_coord*4))/4) + # make it square
    scale_y_continuous(limits = c(min_coord, max_coord), breaks = seq(-ceiling(max_coord*4), ceiling(max_coord*4))/4) +
    annotate("text", x = max_coord, y = max_coord, label = num_pos_pos, size = 3) +
    annotate("text", x = max_coord, y = min_coord, label = num_pos_neg, size = 3) +
    annotate("text", x = min_coord, y = max_coord, label = num_neg_pos, size = 3) +
    annotate("text", x = min_coord, y = min_coord, label = num_neg_neg, size = 3) +
    annotate("text", x = max_coord, y = 0, label = num_pos_not, size = 3) +
    annotate("text", x = min_coord, y = 0, label = num_neg_not, size = 3) +
    annotate("text", x = 0, y = max_coord, label = num_not_pos, size = 3) +
    annotate("text", x = 0, y = min_coord, label = num_not_neg, size = 3) +
    labs(x = sprintf("Delta PSI\n%s", paste0(strsplit(study_1, "_")[[1]], collapse = " ")),
         y = sprintf("%s\nDelta PSI", paste0(strsplit(study_2, "_")[[1]], collapse = " "))) + #,
    # subtitle = sprintf("Rsquared = %s", round(rsq, 4))) +

    geom_abline(intercept = 0, slope = correlation_value, alpha = 0.5, linetype = "F1", color = "dodgerblue2", size = 0.8) +
    theme_bw() +
    theme(axis.text.y   = element_text(size = 8),
          axis.text.x   = element_text(size = 8),
          axis.title.y  = element_text(size = 10),
          axis.title.x  = element_text(size = 10),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1)
    )


  # add the annotations
  if(!slope){
    p <- p +
      annotate("text", x = min_coord + 0.2,
               y = max_coord - 0.2,
               label = substitute(italic(R)^2~"="~r2,
                                  list(r2 = format(rsq, digits = 3))),
               size = 3)
  } else {  # add slope
    p <- p +
      annotate("text", x = min_coord + 0.2,
               y = max_coord - 0.2,
               label = substitute(italic(R)^2~"="~r2,
                                  list(r2 = format(rsq, digits = 3))),
               size = 3) +
      annotate("text", x = min_coord + 0.2,
               y = max_coord - 0.30,
               label = substitute("m ="~slop,
                                  list(slop = format(correlation_value, digits = 3))),
               size = 3)


  }

  return(p)
}
