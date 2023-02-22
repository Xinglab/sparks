###### FUNCTIONS  #####
#' @export
calculate_RBP_KD_correlation_from_mats <- function(kd_mats_df,
                                                   study_sig_df,
                                                   study_1,
                                                   study_2,
                                                   # perm_mats,
                                                   beta_threshold = 0.1,
                                                   num_permutation = 100,
                                                   perm = F){

  # merge each mats by event
  merged_mats_df <- merge(kd_mats_df, study_sig_df, by = "event")

  delta_psi_df_prefilter <- merged_mats_df[, c("beta.x", "beta.y")]
  rownames(delta_psi_df_prefilter) <- merged_mats_df$event
  # delta_psi_df_prefilter <- merged_mats_df[abs(merged_mats_df$beta.x) > 0.1 | merged_mats_df$fdr.x < 0.01, c("beta.x", "beta.y")]

  # delta_psi_df <- merged_mats_df[merged_mats_df$fdr.x < 0.01 | merged_mats_df$fdr.y < 0.01, c("beta.x", "beta.y")]

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
  if(is.na(score)){ score <- 0}  # fill 0 if score is NaN, which is possible by 0 div
  score_abs <- abs(score)


  ## Permutation for significance calculation
  # # generate permutation table
  # random_index <- sample(ncol(perm_mats_filtered), num_permutation)
  # score_distribution <- unlist(lapply(seq(num_permutation), function(idx){
  #   # swap permutation value
  #   kd_mats_df_perm <- kd_mats_df
  #   perm_beta <- perm_mats[, random_index[idx]]
  #
  #   kd_mats_df_perm$beta <- perm_beta
  #
  #   merged_mats_df <- merge(kd_mats_df_perm, study_sig_df, by = "event")
  #
  #   delta_psi_df_prefilter <- merged_mats_df[, c("beta.x", "beta.y")]
  #   colnames(delta_psi_df_prefilter) <- c("S1", "S2")
  #   # keep only complete observations
  #   delta_psi_df <- delta_psi_df_prefilter[rowSums(is.na(delta_psi_df_prefilter)) == 0, ]
  #
  #   perm_score <- calculate_permuted_score_from_psi_df(delta_psi_df)
  #   return(perm_score)
  # }))
  # TEST - original model
  # IF PERM = TRUE
  if (perm){


    score_distribution <- unlist(lapply(seq(num_permutation), function(idx){

      perm_score <- calculate_permuted_score_from_psi_df(psi_df_any_sig, beta_threshold = beta_threshold)
      return(perm_score)
    }))

    # calculate significant percentile threshold for future comparisons
    score_quantiles <- quantile(score_distribution, c(0.025, 0.975), na.rm = T)
    score_quantiles_neg <- score_quantiles[1]
    score_quantiles_pos <- score_quantiles[2]

    # calculate perm pval
    perm_pval_pos <- (sum(score_distribution > score)) / num_permutation
    perm_pval_neg <- (sum(score_distribution < score)) / num_permutation

    # choose the lower pval - as this is bidirectional
    perm_pval <-  min(perm_pval_pos, perm_pval_neg)
    perm_pval_abs <- sum(abs(score_distribution) > score_abs) / num_permutation
    result_df <- data.frame(S1 = study_1,
                            S2 = study_2,
                            correlation = correlation_value,
                            rsquared = rsq,
                            num_common_events = num_complete,
                            num_S1_events = num_S1_sig,
                            num_S2_events = num_S2_sig,
                            num_pos_pos = num_pos_pos,
                            num_neg_neg = num_neg_neg,
                            num_pos_neg = num_pos_neg,
                            num_neg_pos = num_neg_pos,
                            num_not_pos = num_not_pos,
                            num_not_neg = num_not_neg,
                            num_pos_not = num_pos_not,
                            num_neg_not = num_neg_not,
                            num_eff = num_eff,
                            num_total_sig = num_total,
                            score = score,
                            score_abs = score_abs,
                            pval = pval,
                            pval_perm = perm_pval,
                            pval_pos = perm_pval_pos,
                            pval_neg = perm_pval_neg,
                            pval_perm_abs = perm_pval_abs,
                            score_threshold_pos = score_quantiles_pos,
                            score_threshold_neg = score_quantiles_neg)
  } else {


    result_df <- data.frame(S1 = study_1,
                            S2 = study_2,
                            correlation = correlation_value,
                            rsquared = rsq,
                            num_common_events = num_complete,
                            num_S1_events = num_S1_sig,
                            num_S2_events = num_S2_sig,
                            num_pos_pos = num_pos_pos,
                            num_neg_neg = num_neg_neg,
                            num_pos_neg = num_pos_neg,
                            num_neg_pos = num_neg_pos,
                            num_not_pos = num_not_pos,
                            num_not_neg = num_not_neg,
                            num_pos_not = num_pos_not,
                            num_neg_not = num_neg_not,
                            num_eff = num_eff,
                            num_total_sig = num_total,
                            score = score,
                            score_abs = score_abs,
                            pval = pval)
  }

  return(result_df)
}



#' @export
calculate_RBP_KD_correlation_from_mats_fdr <- function(kd_mats_df,
                                                       study_sig_df,
                                                       study_1,
                                                       study_2,
                                                       # perm_mats,
                                                       beta_threshold = 0.1,
                                                       num_permutation = 100,
                                                       perm = F,
                                                       fdr_threshold = 0.05){


  # merge each mats by event
  merged_mats_df <- merge(kd_mats_df, study_sig_df, by = "event")

  # delta_psi_df_prefilter <- merged_mats_df[, c("beta.x", "beta.y")]
  # rownames(delta_psi_df_prefilter) <- merged_mats_df$event
  # delta_psi_df_prefilter <- merged_mats_df[abs(merged_mats_df$beta.x) > 0.1 | merged_mats_df$fdr.x < 0.01, c("beta.x", "beta.y")]

  delta_psi_df_prefilter <- merged_mats_df[merged_mats_df$fdr.x < fdr_threshold | merged_mats_df$fdr.y < fdr_threshold, c("beta.x", "beta.y")]

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


  ## Permutation for significance calculation
  # # generate permutation table
  # random_index <- sample(ncol(perm_mats_filtered), num_permutation)
  # score_distribution <- unlist(lapply(seq(num_permutation), function(idx){
  #   # swap permutation value
  #   kd_mats_df_perm <- kd_mats_df
  #   perm_beta <- perm_mats[, random_index[idx]]
  #
  #   kd_mats_df_perm$beta <- perm_beta
  #
  #   merged_mats_df <- merge(kd_mats_df_perm, study_sig_df, by = "event")
  #
  #   delta_psi_df_prefilter <- merged_mats_df[, c("beta.x", "beta.y")]
  #   colnames(delta_psi_df_prefilter) <- c("S1", "S2")
  #   # keep only complete observations
  #   delta_psi_df <- delta_psi_df_prefilter[rowSums(is.na(delta_psi_df_prefilter)) == 0, ]
  #
  #   perm_score <- calculate_permuted_score_from_psi_df(delta_psi_df)
  #   return(perm_score)
  # }))
  # TEST - original model
  # IF PERM = TRUE
  if (perm){


    score_distribution <- unlist(lapply(seq(num_permutation), function(idx){

      perm_score <- calculate_permuted_score_from_psi_df(psi_df_any_sig, beta_threshold = beta_threshold)
      return(perm_score)
    }))

    # calculate significant percentile threshold for future comparisons
    score_quantiles <- quantile(score_distribution, c(0.025, 0.975), na.rm = T)
    score_quantiles_neg <- score_quantiles[1]
    score_quantiles_pos <- score_quantiles[2]

    # calculate perm pval
    perm_pval_pos <- (sum(score_distribution > score)) / num_permutation
    perm_pval_neg <- (sum(score_distribution < score)) / num_permutation

    # choose the lower pval - as this is bidirectional
    perm_pval <-  min(perm_pval_pos, perm_pval_neg)
    perm_pval_abs <- sum(abs(score_distribution) > score_abs) / num_permutation
    result_df <- data.frame(S1 = study_1,
                            S2 = study_2,
                            correlation = correlation_value,
                            rsquared = rsq,
                            num_common_events = num_complete,
                            num_S1_events = num_S1_sig,
                            num_S2_events = num_S2_sig,
                            num_pos_pos = num_pos_pos,
                            num_neg_neg = num_neg_neg,
                            num_pos_neg = num_pos_neg,
                            num_neg_pos = num_neg_pos,
                            num_not_pos = num_not_pos,
                            num_not_neg = num_not_neg,
                            num_pos_not = num_pos_not,
                            num_neg_not = num_neg_not,
                            num_eff = num_eff,
                            num_total_sig = num_total,
                            score = score,
                            score_abs = score_abs,
                            pval = pval,
                            pval_perm = perm_pval,
                            pval_pos = perm_pval_pos,
                            pval_neg = perm_pval_neg,
                            pval_perm_abs = perm_pval_abs,
                            score_threshold_pos = score_quantiles_pos,
                            score_threshold_neg = score_quantiles_neg)
  } else {


    result_df <- data.frame(S1 = study_1,
                            S2 = study_2,
                            correlation = correlation_value,
                            rsquared = rsq,
                            num_common_events = num_complete,
                            num_S1_events = num_S1_sig,
                            num_S2_events = num_S2_sig,
                            num_pos_pos = num_pos_pos,
                            num_neg_neg = num_neg_neg,
                            num_pos_neg = num_pos_neg,
                            num_neg_pos = num_neg_pos,
                            num_not_pos = num_not_pos,
                            num_not_neg = num_not_neg,
                            num_pos_not = num_pos_not,
                            num_neg_not = num_neg_not,
                            num_eff = num_eff,
                            num_total_sig = num_total,
                            score = score,
                            score_abs = score_abs,
                            pval = pval)
  }

  return(result_df)
}



# permute PSI df and calculate score
#' @export
calculate_permuted_score_from_psi_df <- function(delta_psi_df,
                                                 beta_threshold = 0.1,
                                                 num_permutation = 100,
                                                 perm_threshold = 0.05){
  # subset events with any sig
  # this will elimiante events w/in (|0.1|, |0.1|)
  psi_df_any_sig <- subset(delta_psi_df, abs(S1) > beta_threshold | abs(S2) > beta_threshold)

  ### Permute the dataframe
  psi_df_any_sig_permuted <- psi_df_any_sig
  ## permute in each significant quadrant only
  # extract the quadrant info
  pos_pos_quad <- subset(psi_df_any_sig_permuted, S1 > perm_threshold & S2 > perm_threshold)
  pos_neg_quad <- subset(psi_df_any_sig_permuted, S1 > perm_threshold & S2 < -perm_threshold)
  neg_pos_quad <- subset(psi_df_any_sig_permuted, S1 < -perm_threshold & S2 > perm_threshold)
  neg_neg_quad <- subset(psi_df_any_sig_permuted, S1 < -perm_threshold & S2 < -perm_threshold)

  # shuffle S2
  psi_df_any_sig_permuted[rownames(pos_pos_quad), ]$S2 <- sample(pos_pos_quad$S2)
  psi_df_any_sig_permuted[rownames(pos_neg_quad), ]$S2 <- sample(pos_neg_quad$S2)
  psi_df_any_sig_permuted[rownames(neg_pos_quad), ]$S2 <- sample(neg_pos_quad$S2)
  psi_df_any_sig_permuted[rownames(neg_neg_quad), ]$S2 <- sample(neg_neg_quad$S2)


  # calculate summary data to find beta direction
  num_pos_pos <- dim(subset(psi_df_any_sig_permuted, S1 > beta_threshold & S2 > beta_threshold))[1]
  num_pos_neg <- dim(subset(psi_df_any_sig_permuted, S1 > beta_threshold & S2 < -beta_threshold))[1]
  num_neg_pos <- dim(subset(psi_df_any_sig_permuted, S1 < -beta_threshold & S2 > beta_threshold))[1]
  num_neg_neg <- dim(subset(psi_df_any_sig_permuted, S1 < -beta_threshold & S2 < -beta_threshold))[1]

  # calculate neutral ones for filtering as well
  num_not_pos <- dim(subset(psi_df_any_sig_permuted, S1 < beta_threshold & S1 > -beta_threshold & S2 > beta_threshold))[1]
  num_not_neg <- dim(subset(psi_df_any_sig_permuted, S1 < beta_threshold & S1 > -beta_threshold & S2 < -beta_threshold))[1]
  num_pos_not <- dim(subset(psi_df_any_sig_permuted, S1 > beta_threshold & S2 < beta_threshold & S2 > -beta_threshold))[1]
  num_neg_not <- dim(subset(psi_df_any_sig_permuted, S1 < -beta_threshold & S2 < beta_threshold & S2 > -beta_threshold))[1]


  psi_df_both_sig <- subset(psi_df_any_sig_permuted, abs(S1) > beta_threshold & abs(S2) > beta_threshold)

  # run linear model and get relevant statistics
  tryCatch({
    # lm_result <- lm(psi_df_any_sig_permuted$S2 ~ psi_df_any_sig_permuted$S1 + 0)
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

  # dev test - permutation
  ## permutation for p-value calculation


  # calculate effective number
  num_eff <- (num_pos_pos + num_neg_neg) - (num_neg_pos + num_pos_neg)
  num_total <- (num_pos_pos + num_neg_neg) + (num_neg_pos + num_pos_neg) + (num_not_pos + num_not_neg) + (num_pos_not + num_neg_not)

  # calculate score
  score <- rsq * num_eff / num_total
  score_abs <- abs(score)
  return(score)
}




#' @export
calculate_RBP_KD_concordance_from_mats <- function(kd_mats_df,
                                                   study_sig_df,
                                                   study_1,
                                                   study_2,
                                                   beta_threshold = 0.1){

  # merge each mats by event
  merged_mats_df <- merge(kd_mats_df, study_sig_df, by = "event")

  delta_psi_df_prefilter <- merged_mats_df[, c("beta.x", "beta.y")]
  # delta_psi_df_prefilter <- merged_mats_df[abs(merged_mats_df$beta.x) > 0.1 | merged_mats_df$fdr.x < 0.01, c("beta.x", "beta.y")]

  # delta_psi_df <- merged_mats_df[merged_mats_df$fdr.x < 0.01 | merged_mats_df$fdr.y < 0.01, c("beta.x", "beta.y")]

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

  # calculate summary data to find beta direction
  num_pos_pos <- dim(subset(delta_psi_df, S1 > beta_threshold & S2 > beta_threshold))[1]
  num_pos_neg <- dim(subset(delta_psi_df, S1 > beta_threshold & S2 < -beta_threshold))[1]
  num_neg_pos <- dim(subset(delta_psi_df, S1 < -beta_threshold & S2 > beta_threshold))[1]
  num_neg_neg <- dim(subset(delta_psi_df, S1 < -beta_threshold & S2 < -beta_threshold))[1]

  # calculate neutral ones for filtering as well
  num_not_pos <- dim(subset(delta_psi_df, S1 < beta_threshold & S1 > -beta_threshold & S2 > beta_threshold))[1]
  num_not_neg <- dim(subset(delta_psi_df, S1 < beta_threshold & S1 > -beta_threshold & S2 < -beta_threshold))[1]
  num_pos_not <- dim(subset(delta_psi_df, S1 > beta_threshold & S2 < beta_threshold & S2 > -beta_threshold))[1]
  num_neg_not <- dim(subset(delta_psi_df, S1 < -beta_threshold & S2 < beta_threshold & S2 > -beta_threshold))[1]

  # calculate effective number
  num_eff <- (num_pos_pos + num_neg_neg) - (num_neg_pos + num_pos_neg)
  num_total <- (num_pos_pos + num_neg_neg) + (num_neg_pos + num_pos_neg) + (num_not_pos + num_not_neg) + (num_pos_not + num_neg_not)

  # calculate score
  score <- num_eff / num_total
  if(is.na(score)) {score <- 0}
  score_abs <- abs(score)

  # calculate p-value
  sign_df <- ifelse(abs(delta_psi_df) < beta_threshold,
                    0, ifelse(delta_psi_df > 0, 1, -1))

  tryCatch({
    p_value <- cor.test(sign_df[, 1], sign_df[, 2], method = "kendall")$p.value
  },
  error = function(e){
    p_value <<- 1
  })
  result_df <- data.frame(S1 = study_1,
                          S2 = study_2,
                          num_common_events = num_complete,
                          num_S1_events = num_S1_sig,
                          num_S2_events = num_S2_sig,
                          num_pos_pos = num_pos_pos,
                          num_neg_neg = num_neg_neg,
                          num_pos_neg = num_pos_neg,
                          num_neg_pos = num_neg_pos,
                          num_not_pos = num_not_pos,
                          num_not_neg = num_not_neg,
                          num_pos_not = num_pos_not,
                          num_neg_not = num_neg_not,
                          num_eff = num_eff,
                          num_total_sig = num_total,
                          score = score,
                          score_abs = score_abs,
                          pval = p_value)

  return(result_df)
}

#' @export
calculate_RBP_KD_concordance_from_mats_fdr <- function(kd_mats_df,
                                                       study_sig_df,
                                                       study_1,
                                                       study_2,
                                                       beta_threshold = 0.1,
                                                       fdr_threshold = 0.05){

  # merge each mats by event
  merged_mats_df <- merge(kd_mats_df, study_sig_df, by = "event")

  # delta_psi_df_prefilter <- merged_mats_df[, c("beta.x", "beta.y")]
  # delta_psi_df_prefilter <- merged_mats_df[abs(merged_mats_df$beta.x) > 0.1 | merged_mats_df$fdr.x < 0.01, c("beta.x", "beta.y")]

  delta_psi_df_prefilter <- merged_mats_df[merged_mats_df$fdr.x < fdr_threshold | merged_mats_df$fdr.y < fdr_threshold, c("beta.x", "beta.y")]

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

  # calculate summary data to find beta direction
  num_pos_pos <- dim(subset(delta_psi_df, S1 > beta_threshold & S2 > beta_threshold))[1]
  num_pos_neg <- dim(subset(delta_psi_df, S1 > beta_threshold & S2 < -beta_threshold))[1]
  num_neg_pos <- dim(subset(delta_psi_df, S1 < -beta_threshold & S2 > beta_threshold))[1]
  num_neg_neg <- dim(subset(delta_psi_df, S1 < -beta_threshold & S2 < -beta_threshold))[1]

  # calculate neutral ones for filtering as well
  num_not_pos <- dim(subset(delta_psi_df, S1 < beta_threshold & S1 > -beta_threshold & S2 > beta_threshold))[1]
  num_not_neg <- dim(subset(delta_psi_df, S1 < beta_threshold & S1 > -beta_threshold & S2 < -beta_threshold))[1]
  num_pos_not <- dim(subset(delta_psi_df, S1 > beta_threshold & S2 < beta_threshold & S2 > -beta_threshold))[1]
  num_neg_not <- dim(subset(delta_psi_df, S1 < -beta_threshold & S2 < beta_threshold & S2 > -beta_threshold))[1]

  # calculate effective number
  num_eff <- (num_pos_pos + num_neg_neg) - (num_neg_pos + num_pos_neg)
  num_total <- (num_pos_pos + num_neg_neg) + (num_neg_pos + num_pos_neg) + (num_not_pos + num_not_neg) + (num_pos_not + num_neg_not)

  # calculate score
  score <- num_eff / num_total
  score_abs <- abs(score)

  result_df <- data.frame(S1 = study_1,
                          S2 = study_2,
                          num_common_events = num_complete,
                          num_S1_events = num_S1_sig,
                          num_S2_events = num_S2_sig,
                          num_pos_pos = num_pos_pos,
                          num_neg_neg = num_neg_neg,
                          num_pos_neg = num_pos_neg,
                          num_neg_pos = num_neg_pos,
                          num_not_pos = num_not_pos,
                          num_not_neg = num_not_neg,
                          num_pos_not = num_pos_not,
                          num_neg_not = num_neg_not,
                          num_eff = num_eff,
                          num_total_sig = num_total,
                          score = score,
                          score_abs = score_abs)
  # sig events in each RBP
  # sig events common
  # corr
  # corr stat
  return(result_df)
}


#' @export
GSEA.EnrichmentScore <-
  function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {

    tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <- N - Nh
    if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector^alpha)
    sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
    norm.tag <- 1/sum.correl.tag
    norm.no.tag <- 1/Nm
    RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      # ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      # ES <- min.ES
      ES <- signif(min.ES, digits = 5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
  }


#' @export
calculate_GSEA_score <- function(input_signature_list, test_pos_events, test_neg_events){
  pos_overlap_list <- test_pos_events[test_pos_events %in% input_signature_list$event]
  if (length(pos_overlap_list) > 0){
    pos_result <- GSEA.EnrichmentScore(input_signature_list$event, pos_overlap_list, weighted.score.type = 0)
    pos_stat <- pos_result$ES
  } else {
    pos_stat <- 0
  }

  neg_overlap_list <- test_neg_events[test_neg_events %in% input_signature_list$event]
  if (length(neg_overlap_list) > 0){
    neg_result <- GSEA.EnrichmentScore(input_signature_list$event, neg_overlap_list, weighted.score.type = 0)
    neg_stat <- neg_result$ES
  } else {
    neg_stat <- 0
  }

  result_df <- data.frame(pos_score = pos_stat,
                          neg_score = neg_stat)
  return(result_df)
}


#' @export
calculate_pulled_binomial_pval <- function(psi_1,
                                           psi_2,
                                           total_count_1,
                                           total_count_2){
  prop_matrix <- t(matrix(c(total_count_1 * psi_1, total_count_1 * (1 - psi_1), total_count_2 * psi_2, total_count_2 * (1 - psi_2)), nrow = 2))

  prop_test_result <- prop.test(prop_matrix, alternative = "two.sided")
  prop_pval <- prop_test_result$p.value
  return(prop_pval)
}


# pseudo binomial pval
#' @export
calculate_pulled_binomial_pval_for_mats <- function(input_mats_df) {
  pval_list <- unlist(apply(input_mats_df, 1, function(x){

    # calcualte PSI
    psi_1 <- x['pulled_psi_1']
    psi_2 <- x['pulled_psi_2']

    count_array <- as.numeric(strsplit(x['count_values'], ",")[[1]])
    total_count_1 <- sum(count_array[1:2])
    total_count_2 <- sum(count_array[3:4])

    p <- calculate_pulled_binomial_pval(as.numeric(psi_1), as.numeric(psi_2), as.numeric(total_count_1), as.numeric(total_count_2))
    return(p)
  }))

  input_mats_df$pulled_pval <- pval_list
  return(input_mats_df)
}

