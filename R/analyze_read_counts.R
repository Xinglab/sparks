###### FUNCTIONS #####
# categorize

#' @export
categorize_mats_df <- function(input_mats, beta_threshold = 0.1, fdr_threshold = 0.01){
  input_mats$beta_sig <- ifelse(abs(input_mats$pulled_delta_psi) > beta_threshold, "sig", "notsig")
  input_mats$fdr_sig <- ifelse(input_mats$fdr < fdr_threshold, "sig", "notsig")
  input_mats$sig <- ifelse(abs(input_mats$pulled_delta_psi) > beta_threshold & input_mats$fdr < fdr_threshold, "sig", "notsig")

  return(input_mats)
}


#' @export
calculate_log2_count_mode <- function(input_mats){
  count_dens <- density(log2(input_mats$avg_count + 1))

  # calculate mode
  count_mode <- count_dens$x[which.max(count_dens$y)]
  # count_mode <- median(count_dens$x)
  return(count_mode)
}


#' @export
calculate_event_stats_for_mats <- function(study_mats, study){
  # categorize by beta and fdr
  cat_study_mats <- categorize_mats_df(study_mats)

  # subset each
  nodiff_entries <- subset(cat_study_mats, beta_sig == "notsig" & fdr_sig == "notsig")

  beta_sig_only_entries <- subset(cat_study_mats, beta_sig == "sig" & fdr_sig == "notsig")

  fdr_sig_only_entries <- subset(cat_study_mats, beta_sig == "notsig" & fdr_sig == "sig")

  sig_entries <- subset(cat_study_mats, beta_sig == "sig" & fdr_sig == "sig")

  # calculate mode
  nodiff_mode <- calculate_log2_count_mode(nodiff_entries)
  beta_sig_mode <- calculate_log2_count_mode(beta_sig_only_entries)
  dummy <- tryCatch({
    fdr_sig_mode <- calculate_log2_count_mode(fdr_sig_only_entries)},
    error=function(e){fdr_sig_mode <<- 0})
  dummy <- tryCatch({
    sig_mode <- calculate_log2_count_mode(sig_entries)},
    error=function(e){sig_mode <<- 0})

  # calculate count dist diff
  tryCatch({
    pval <- t.test(log2(beta_sig_only_entries$avg_count + 1),
                   log2(sig_entries$avg_count + 1))$p.value
  },
  error = function(e){pval <<- 1})

  # aggregate data and return

  result_df <- data.frame(study = study,
                          pval = pval,
                          nodiff_mode = nodiff_mode,
                          beta_sig_mode = beta_sig_mode,
                          fdr_sig_mode = fdr_sig_mode,
                          sig_mode = sig_mode,
                          nodiff_count = dim(nodiff_entries)[1],
                          beta_sig_count = dim(beta_sig_only_entries)[1],
                          fdr_sig_count = dim(fdr_sig_only_entries)[1],
                          sig_count = dim(sig_entries)[1])
  return(result_df)
}


