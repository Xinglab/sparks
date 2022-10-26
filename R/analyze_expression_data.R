# DIFFERENTIAL EXPRESSION


# TODO - add documentation here
perform_DE_analysis <- function(object, group_one, group_two){

  # extract expression data
  expression_df <- object@exp_df

  diff_exp_result_df <- do.call(rbind, apply(expression_df, 1, function(exp){
    exp_one <- exp[group_one]
    exp_two <- exp[group_two]

    # calculate pval
    sig_test <- t.test(exp_one, exp_two, alternative = "two.sided")
    pval <- sig_test$p.value

    # caluclate beta
    delta_tpm <- sig_test$estimate[1] - sig_test$estimate[2]
    log2fc <- log2(sig_test$estimate[1]/sig_test$estimate[2])

    # gather results
    result_df <- data.frame(beta = delta_tpm,
                            log2fc = log2fc,
                            mean_exp_group_one = mean(exp_one, na.rm = T),
                            mean_exp_group_two = mean(exp_two, na.rm = T),
                            pval = pval)
    return(result_df)

  }))
  # drop NA
  diff_exp_result_df_filtered <- diff_exp_result_df[rowSums(is.na(diff_exp_result_df)) == 0, ]

  # calculate FDR
  diff_exp_result_df_filtered$qval <- p.adjust(diff_exp_result_df_filtered$pval, method = "BH")

  # update the slot
  object@diff_exp_result <- diff_exp_result_df_filtered
  return(object)
}


query_expression_data_for_gene_set <- function(object, gene_list){
  expression_df <- object@exp_df

  # extract shortened gene name
  # trimmed_gene_name <- unlist(lapply(rownames(expression_df), function(x) strsplit(x, ":")[[1]][2]))
  # rownames(expression_df) <- trimmed_gene_name
  # expression_df$gene <- trimmed_gene_name
  # match with the given gene list
  # target_exp_df <- expression_df[gene_list,, drop = F]


  trimmed_gene_name <- extract_gene_symbols(rownames(expression_df))

  target_exp_df <- expression_df[trimmed_gene_name %in% gene_list, ]
  # add the gene for melting for downstream plotting
  target_exp_df$gene <- extract_gene_symbols(rownames(target_exp_df))
  rownames(target_exp_df) <- NULL
  return(target_exp_df)
}


# TODO - add documentation
remove_expression_duplicates <- function(object){
  trimmed_gene_name <- unlist(lapply(rownames(object@exp_df), function(x) strsplit(x, ":")[[1]][2]))


  not_duplicates <- object@exp_df[!duplicated(trimmed_gene_name) & !duplicated(trimmed_gene_name, fromLast = T), ]

  duplicates <- object@exp_df[duplicated(trimmed_gene_name)| duplicated(trimmed_gene_name, fromLast = T), ]
  duplicates$gene <- unlist(lapply(rownames(duplicates), function(x) strsplit(x, ":")[[1]][2]))

  deduped_exp_df <- do.call(rbind, lapply(unique(duplicates$gene), function(Gene){
    gene_dup <- subset(duplicates, gene == Gene)
    gene_dup$gene <- NULL

    # collapse by sum
    gene_collapsed <- as.data.frame(t(colSums(gene_dup)))
    rownames(gene_collapsed) <- rownames(gene_dup)[1]
    # gene_collapsed$gene <- Gene

    return(gene_collapsed)
  }))

  clean_exp_df <- rbind(not_duplicates, deduped_exp_df)

  # update back to the object
  object@exp_df <- clean_exp_df
  return(object)
}


# TODO - add documentation
generate_expression_boxplot_for_DE <- function(object, group_one, group_two, gene_list,
                                               group_one_label = "Group One",
                                               group_two_label = "Group Two",
                                               group_one_color = "indianred3",
                                               group_two_color = "pink1",
                                               sample_description = "Engineered Cell line - MYC"){
  target_exp_df <- query_expression_data_for_gene_set(object, gene_list)
  plot_df <- reshape2::melt(target_exp_df)
  # sample_levels <- c( "ICA1.on","ICA1.off","ICA2.on","ICA2.off", "ICA3.on","ICA3.off"
  sample_levels <- c(group_one, group_two)
  plot_df$variable <- factor(plot_df$variable, levels = sample_levels)

  # set labels
  label_vector <- c(rep(group_one_label, length(group_one)), rep(group_two_label, length(group_two)))
  names(label_vector) <- sample_levels
  plot_df$label <- label_vector[plot_df$variable]


  p <- ggplot(plot_df, aes( x = label, y = value, fill = label)) +
    geom_boxplot() + geom_jitter(width = 0.1) +
    labs(y = sprintf("Expression (TPM)")) +
    theme(legend.position = "none") +
    scale_fill_manual(values = setNames(c(group_one_color, group_two_color), c(group_one_label, group_two_label))) +
    facet_wrap(~gene, scales = "free_y") +
    labs(x = sample_description)

  return(p)
}





