#'@export
generate_enrichment_barplot <- function(result_df,
                                        bar_color = "grey75",
                                        sig_threshold = 0.01,
                                        rank_label = T,
                                        sig_test_method = "bonferroni",
                                        num_plot = NA,
                                        select_genes = c(),
                                        manual_colors = list(),
                                        text_scale_factor = 1,
                                        select_gene_marker = F){
  # calculate padj if it isn't there
  if (!("padj" %in% colnames(result_df))){
    # print("padj missing, calculating again")
    result_df <- calculate_SPARKS_padj(result_df)
  }

  # calculate rank
  result_df$plot_rank <- rank(-result_df$gsea_score)

  # define genes for downstream processing
  result_df$rbp_s1 <- unlist(lapply(result_df$S1, function(x) strsplit(x, "_")[[1]][2]))

  # subset the significant ones
  sig_subset_df <- subset(result_df, padj < sig_threshold)

  # sort by score
  sig_subset_df$S1 <- factor(sig_subset_df$S1,
                             levels = rev((sig_subset_df %>% arrange(gsea_score))$S1))

  x_lab <- "Significantly enriched experiments"

  if (length(select_genes) == 0){

    # subset if top N is given and change plot annotation
    if (!is.na(num_plot)){
      # change x lab
      x_lab <- sprintf("Top %s significantly enriched experiments\nTotal # significant experiment = %s", min(num_plot, dim(sig_subset_df)[1]), dim(sig_subset_df)[1])

      # subset by top N
      sig_subset_df <- sig_subset_df %>% top_n(abs(gsea_score), n = num_plot)
    }

    # generate plot
    p <-
      ggplot(sig_subset_df) +
      geom_bar(aes(x = S1,
                   y = gsea_score),
               stat = 'identity',
               fill = bar_color,
               color = "black") +
      geom_text(aes(x = S1,
                    label = lapply(S1, function(x) gsub(x, pattern = "_", replacement = " ")),
                    y = ifelse(gsea_score > 0, -0.05, 0.05),
                    hjust = ifelse(gsea_score > 0, 1, 0),
                    color = ifelse((rbp_s1 %in% names(manual_colors)) & (length(names(manual_colors)) > 0),
                                   manual_colors[rbp_s1],
                                   "black")),
                                   angle = 90,
                size = 3 * text_scale_factor)+
      theme(axis.text.y   = element_text(size = 8 * text_scale_factor),
            axis.text.x   = element_text(size = 8 * text_scale_factor),
            axis.title.y  = element_text(size = 10 * text_scale_factor),
            axis.title.x  = element_text(size = 10 * text_scale_factor),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # axis.line = element_line(colour = "black"),
            panel.border = element_rect(color = "black",
                                        fill = NA,
                                        size = 1),
            legend.position = "bottom"
      ) +
      scale_y_continuous(limits = c(-ceiling(max(abs(sig_subset_df$gsea_score) + 0.1) * 4),
                                    ceiling(max(abs(sig_subset_df$gsea_score) + 0.1) * 4)) / 4,
                         breaks = seq(-ceiling(max(abs(sig_subset_df$gsea_score) + 0.1) * 4),
                                      ceiling(max(abs(sig_subset_df$gsea_score) + 0.1) * 4)) / 4) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(x = x_lab,
           y = "SPARKS Enrichment Score") +
    scale_color_identity()


    # add rank label on top
    if(rank_label){
      p <-
        p + geom_text(aes(label = plot_rank,
                          x = S1,
                          y = ifelse(gsea_score > 0, gsea_score + 0.1, gsea_score - 0.1)),
                      size = 2.25 * text_scale_factor,
                      angle = 90,
                      vjust = 0.5)
    }

    return(p)
  } else {  # if there are genes of interest
    # choose ranks in the select rbps
    select_gene_result_df <- subset(result_df, rbp_s1 %in% select_genes)

    # subset if top N is given and change plot annotation
    if (!is.na(num_plot)){
      # change x lab
      x_lab <- sprintf("Select experiments\nTotal # significant experiment = %s", dim(sig_subset_df)[1])

      # subset by top N
      sig_subset_df <- sig_subset_df %>% top_n(abs(gsea_score), n = num_plot)

    }

    # merge the result and take unique to remove duplicates
    merged_subset_df <- unique(rbind(select_gene_result_df, sig_subset_df))

    # sort by score
    merged_subset_df$S1 <- factor(merged_subset_df$S1,
                                  levels = merged_subset_df$S1[order(-merged_subset_df$gsea_score)])

    # denote p-val sig and rank sig
    merged_subset_df$rank_sig <- ifelse(merged_subset_df$plot_rank > 0.05 * max(result_df$plot_rank) &
                                          merged_subset_df$plot_rank < 0.95 * max(result_df$plot_rank), "notsig", "sig")
    merged_subset_df$pval_sig <- ifelse(merged_subset_df$padj > sig_threshold, "notsig", "sig")

    p <-
      ggplot(merged_subset_df) +
      geom_bar(aes(x = S1,
                   y = gsea_score),
               stat = 'identity',
               fill = bar_color,
               color = "black") +
      geom_text(aes(x = S1,
                    label = lapply(S1, function(x) gsub(x, pattern = "_", replacement = " ")),
                    y = ifelse(gsea_score > 0, -0.05, 0.05),
                    hjust = ifelse(gsea_score > 0, 1, 0),
                    color = ifelse((rbp_s1 %in% names(manual_colors)) & (length(names(manual_colors)) > 0),
                                   manual_colors[rbp_s1],
                                   "black")),
                                   angle = 90,
                size = 3 * text_scale_factor) +
      theme(axis.text.y   = element_text(size = 8 * text_scale_factor),
            axis.text.x   = element_text(size = 8 * text_scale_factor),
            axis.title.y  = element_text(size = 10 * text_scale_factor),
            axis.title.x  = element_text(size = 10 * text_scale_factor),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # axis.line = element_line(colour = "black"),
            panel.border = element_rect(color = "black",
                                        fill = NA,
                                        size = 1),
            legend.position = "bottom"
      ) +
      scale_y_continuous(limits = c(-ceiling(max(abs(merged_subset_df$gsea_score) + 0.1) * 4),
                                    ceiling(max(abs(merged_subset_df$gsea_score) + 0.1) * 4)) / 4,
                         breaks = seq(-ceiling(max(abs(merged_subset_df$gsea_score) + 0.1) * 4),
                                      ceiling(max(abs(merged_subset_df$gsea_score) + 0.1) * 4)) / 4) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(x = x_lab,
           y = "SPARKS Enrichment Score") +
      scale_color_identity()
    # add sig info
    p <- p + geom_text(data = merged_subset_df,
                       aes(label = ifelse(pval_sig == "sig", "*", ""),
                           x = S1),
                       y = ceiling(max(abs(merged_subset_df$gsea_score) + 0.1) * 4) / 4 + 0.05,
                       size = 4 * text_scale_factor,
                       vjust = 0.5,
                       hjust = 0,
                       color = "dodgerblue1") +
      geom_text(data = merged_subset_df,
                aes(label = ifelse(rank_sig == "sig", "#", ""),
                    x = S1),
                y = ceiling(max(abs(merged_subset_df$gsea_score) + 0.1) * 4) / 4 + 0.05,
                size = 3 * text_scale_factor,
                vjust = 0,
                hjust = 1,
                color = "indianred2")

    if(rank_label){
      p <-
        p + geom_text(aes(label = plot_rank,
                          x = S1,
                          y = ifelse(gsea_score > 0, gsea_score + 0.1, gsea_score - 0.1)),
                      size = 2.25 * text_scale_factor,
                      angle = 90,
                      vjust = 0.5)
    }

    if(select_gene_marker){
      p <- p +
        labs(subtitle = c("# : within top/bottom 5%\n* : FWER < 0.01"))
    }

    return(p)
  }
}
