# TODO - add DOC
calculate_RBP_KD_correlation_data <- function(beta_df, study_1, study_2){
  delta_psi_df_prefilter <- beta_df[, c(study_1, study_2), drop = F]
  colnames(delta_psi_df_prefilter) <- c("S1", "S2")

  # note number of complete observations
  num_obs <- sum(rowSums(!is.na(delta_psi_df_prefilter)) == 2 )
  num_S1 <- sum(!is.na(delta_psi_df_prefilter$S1))
  num_S2 <- sum(!is.na(delta_psi_df_prefilter$S2))


  # only select complete observations
  delta_psi_df <- delta_psi_df_prefilter[rowSums(is.na(delta_psi_df_prefilter)) == 0, ]

  # calculate summary data to find increased or decreased or discordant
  num_positive_positive <- dim(subset(delta_psi_df, S1 > 0 & S2 > 0))[1]
  num_positive_negative <- dim(subset(delta_psi_df, S1 > 0 & S2 < 0))[1]
  num_negative_positive <- dim(subset(delta_psi_df, S1 < 0 & S2 > 0))[1]
  num_negative_negative <- dim(subset(delta_psi_df, S1 < 0 & S2 < 0))[1]

  # run linear model and get relavent statistics
  tryCatch({
    lm_result <- lm(delta_psi_df$S2 ~ delta_psi_df$S1 + 0)
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
  # store data


  # calculate effective number - this will be essentially concordant vs. discordant
  # - however, all discordant means strong negative correlation, so any direction should work
  # - the score derived from this value will be taken as absolute, so direction is fine for prioritization
  effective_num <- (num_positive_positive + num_negative_negative) - (num_positive_negative + num_negative_positive)

  # concordance index
  concordance <- effective_num / num_obs
  if (abs(effective_num) < 10){
    concordance <- concordance * 0.01
  }

  # calculate "jaccard" index score
  effective_jaccard <- effective_num / (num_S1 + num_S2 - num_obs)

  # calculate score
  score <- effective_jaccard * rsq
  score_abs <- abs(score)



  result_df <- data.frame(S1 = study_1,
                          S2 = study_2,
                          correlation = correlation_value,
                          rsquared = rsq,
                          num_common_events = num_obs,
                          num_S1_events = num_S1,
                          num_S2_events = num_S2,
                          num_pos_pos = num_positive_positive,
                          num_pos_neg = num_positive_negative,
                          num_neg_pos = num_negative_positive,
                          num_neg_neg = num_negative_negative,
                          num_eff = effective_num,
                          concordance = concordance,
                          jaccard = effective_jaccard,
                          score = score,
                          score_abs = score_abs,
                          pval = pval)
  # sig events in each RBP
  # sig events common
  # corr
  # corr stat
  return(result_df)
}





# TODO - add DOC
generate_PSI_correlation_plot <- function(beta_df, study_one, study_two, label_y = "bottom"){
  # extract relevant columns
  select_beta_df <- beta_df[, c(study_one, study_two)]

  # select only complete observations
  select_beta_df_filtered <- select_beta_df[rowSums(is.na(select_beta_df)) == 0, ]
  colnames(select_beta_df_filtered) <- c(study_one, study_two)


  # generate string for linear model stats
  lm_formula <-  y ~ x + 0

  # generate plot
  library(ggpmisc)

  p <- ggplot(select_beta_df_filtered, aes_string(x = study_one, y = study_two)) +
    geom_point() +
    stat_smooth(geom = "line",
                method = "lm",
                se = FALSE,
                color = "indianred1",
                formula = lm_formula,
                linetype = "F1",
                alpha = 0.9, size = 0.5,
                fullrange = T) + # fix the intercept at 0
    xlim(min(select_beta_df_filtered, na.rm = T), max(select_beta_df_filtered, na.rm = T)) +
    ylim(min(select_beta_df_filtered, na.rm = T), max(select_beta_df_filtered, na.rm = T)) +
    stat_poly_eq(aes(label = paste(..eq.label.., "x~~", ..rr.label..,sep = "~")),
                 label.x = "right", label.y = label_y,
                 formula = lm_formula, parse = T, size = 4, color = "indianred3")
  return(p)
}


# TODO - add DOC
generate_correlation_plot_study_vs_knockdown <- function(beta_df, study, KD_experiment, study_beta,
                                                         study_label_on, study_label_off, label_y = "bottom"){
  # gather relevant data
  # betas <- beta_df[, KD_experiment, drop = F]
  # delta_psi_df <- cbind(study_beta, betas)
  # colnames(delta_psi_df) <- c(study, KD_experiment)

  # merge using data table
  beta_dt_study <- as.data.table(study_beta)
  beta_dt_study$event <- rownames(study_beta)

  beta_dt_KD <- as.data.table(beta_df[, KD_experiment])
  beta_dt_KD$event <- rownames(beta_df)

  # merge study and kd expierment beta
  delta_psi_dt <- data.table::merge.data.table(beta_dt_study, beta_dt_KD, by = "event", all = T)

  # mutate into dataframe with rownames in the event
  delta_psi_df <- as.data.frame(delta_psi_dt)
  rownames(delta_psi_df) <- delta_psi_df$event
  delta_psi_df$event <- NULL
  colnames(delta_psi_df) <- c(study, KD_experiment)

  # extract relevant information for KD experiment
  knockdown_cellline <- strsplit(KD_experiment, "_")[[1]][1]
  knockdown_study <- strsplit(KD_experiment, "_")[[1]][2]
  knockdown_strategy <- strsplit(KD_experiment, "_")[[1]][3]



  p <- generate_PSI_correlation_plot(delta_psi_df, study, KD_experiment, label_y = label_y)

  p_annotated <- p +
    labs(y = bquote(paste(.(knockdown_cellline), " ", .(knockdown_study), " ", .(knockdown_strategy), ": ", psi [Knockdown] - psi [Control])),
         # x = expression(paste(study, ": ",  psi [MycOn] - psi [MycOff]))) #, "-", expression(paste(psi [MycOff]))))# " - Myc"[Off])))#, paste("Delta PSI in ICA"))))
         x = bquote(paste(.(study), ": ",  psi [ .(study_label_on)] - psi [.(study_label_off)]))) #, "-", expression(paste(psi [MycOff]))))# " - Myc"[Off])))#, paste("Delta PSI in ICA"))))

  return(p_annotated)
}


##### SUMMARY PLOT FUNCTION #####
generate_knockdown_summary_plot <- function(cor_plot_df, encode_KD_class_df){
  cor_plot_df$tier <- encode_KD_class_df[cor_plot_df$S2, ]$class


  # generate bar plot for cor test
  cor_plot_df$S2 <- factor(cor_plot_df$S2,
                           # levels = rev(cor_plot_df$S2[order(cor_plot_df$qval)]))
                           levels = rev(cor_plot_df$S2[order(cor_plot_df$tier, -cor_plot_df$score_abs)]))

  p_num_ica <- ggplot(cor_plot_df, aes(x = S2, y = num_common_events, fill = tier)) +
    geom_bar(stat ="identity", color = "ivory1") +
    scale_fill_manual(values = c("T1" = "dodgerblue2",
                                 "T2" = "lightsteelblue3",
                                 "T3" = "lightsteelblue2")) +
    labs(fill = "ENCODE RBP KD Classification")


  ica_cor_plot_df <- cor_plot_df[, c("S2", "correlation", "rsquared", "qval", "score_abs")]

  # manipulate data to match
  ica_cor_plot_melt <- reshape2::melt(ica_cor_plot_df)
  ica_cor_plot_melt$value_round <- round(ica_cor_plot_melt$value, 2)
  target_parsed <- unlist(lapply(ica_cor_plot_melt$S2, function(x) paste0(strsplit(as.character(x), "_")[[1]][1:3], collapse = " ")))
  ica_cor_plot_melt$target_parsed <- factor(target_parsed,
                                            levels = rev(target_parsed[order(cor_plot_df$score_abs)]))
  ica_cor_plot_melt$variable <- factor(ica_cor_plot_melt$variable,
                                       levels = rev(c("correlation", "rsquared", "qval", "score_abs")))

  p_cor_ica_beta <- ggplot(subset(ica_cor_plot_melt, variable == "correlation"), aes(x = S2, y = variable, fill = value)) +
    geom_tile(width = 0.9, height = 0.9) +
    scale_fill_distiller(palette = "PRGn", direction = -1) +
    geom_text(aes(label = value_round), size = 3) +
    scale_y_discrete(labels = c("correlation" = "Linear Coefficient", "rsquared" = "R2"), expand = c(0,0)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x,
                                                            width = 15))
  p_cor_ica_r2 <- ggplot(subset(ica_cor_plot_melt, variable == "rsquared"), aes(x = S2, y = variable, fill = value)) +
    geom_tile(width = 0.9, height = 0.9) +
    scale_fill_distiller(palette = "Oranges", direction = 1, limits = c(0, 1)) +
    geom_text(aes(label = value_round), size = 3) +
    scale_y_discrete(labels = c("correlation" = "Linear\nCoefficient", "rsquared" = "R2"), expand = c(0,0)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x,
                                                            width = 15))

  p_cor_ica_fdr <- ggplot(subset(ica_cor_plot_melt, variable == "qval"), aes(x = S2, y = variable, fill = -log10(value))) +
    geom_tile(width = 0.9, height = 0.9) +
    scale_fill_distiller(palette = "PuBu", direction = 1, limits = c(0, max(-log10(subset(ica_cor_plot_melt, variable == "qval")$value)))) +
    geom_text(aes(label = round(-log10(value), digits = 2)), size = 3) +
    scale_y_discrete(labels = c("correlation" = "Linear\nCoefficient", "rsquared" = "R2", "qval" = "-log10(FDR)"), expand = c(0,0)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x,
                                                            width = 15))

  p_cor_ica_score <- ggplot(subset(ica_cor_plot_melt, variable == "score_abs"), aes(x = S2, y = variable, fill = value * 100)) +
    geom_tile(width = 0.9, height = 0.9) +
    scale_fill_distiller(palette = "RdPu", direction = 1, limits = c(0, max(subset(ica_cor_plot_melt, variable == "score_abs")$value * 100))) +
    geom_text(aes(label = round(value* 100, digits = 2) ), size = 3) +
    scale_y_discrete(labels = c("correlation" = "Linear\nCoefficient", "rsquared" = "R2", "qval" = "-log10(FDR)", "score_abs" = "Score"), expand = c(0,0)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x,
                                                            width = 15))





  combined_plot <- cowplot::plot_grid(
    p_num_ica  +
      coord_flip() +
      labs(y = "# Myc associated AS Events\nalso associated in KD model",
           x = "Knockdown Model") +
      theme(legend.position = "bottom",
            legend.direction = "vertical"),
    p_cor_ica_beta +
      theme(axis.title = element_blank(),
            axis.text.y = element_blank()) +
      theme(legend.position = "none") +
      coord_flip(),
    p_cor_ica_r2 +
      theme(legend.position = "none") +
      coord_flip() +
      theme(axis.text.y = element_blank(),
            axis.title = element_blank(),
            axis.ticks.y = element_blank()),
    p_cor_ica_fdr +
      theme(legend.position = "none") +
      coord_flip() +
      theme(axis.text.y = element_blank(),
            axis.title = element_blank(),
            axis.ticks.y = element_blank()),
    p_cor_ica_score +
      theme(legend.position = "none") +
      coord_flip() +
      theme(axis.text.y = element_blank(),
            axis.title = element_blank(),
            axis.ticks.y = element_blank()),
    nrow=1,
    align = "h",
    axis = "tb", rel_widths = c(3, 0.5, 0.5, 0.5, 0.5))
  return(combined_plot)
}
