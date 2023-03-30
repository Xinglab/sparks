library(enrichR)

num_terms = 5

# get the raw df
raw_benchmark_df <- extract_benchmark_df(sig_weighted_info_se$raw_benchmark_list[[6]])

# extract genes that satisfies both rank and sig or neither
rep_count <- table(subset(raw_benchmark_df, gsea_rank <= 10 & gsea_padj < 0.01)$rbp_s2)
nonrep_count <- table(subset(raw_benchmark_df, gsea_rank >= 10 & gsea_padj > 0.01)$rbp_s2)

rep_genes <- names(rep_count)[rep_count == 2]
nonrep_genes <- names(nonrep_count)[nonrep_count == 2]

# extract experiment list for query
rep_exp_list <- subset(raw_benchmark_df, rbp_s2 %in% rep_genes)$S2
nonrep_exp_list <- subset(raw_benchmark_df, rbp_s2 %in% nonrep_genes)$S2

# read in the library
kd_library_se <- readRDS("/Users/harryyang/transfer/Clean_minfiltered.SE.library.rds")

# extract the sig count
rep_sig_counts <- do.call(rbind, lapply(rep_exp_list, function(input_exp){
  sig_event_list <- unlist(extract_GSEA_significant_events(kd_library_se[[input_exp]]))
  result_df <- data.frame(exp = input_exp,
                          sig_count = length(sig_event_list),
                          sig = "rep")
  return(result_df)
}))
nonrep_sig_counts <- do.call(rbind, lapply(nonrep_exp_list, function(input_exp){
  sig_event_list <- unlist(extract_GSEA_significant_events(kd_library_se[[input_exp]]))
  result_df <- data.frame(exp = input_exp,
                          sig_count = length(sig_event_list),
                          sig = "nonrep")
  return(result_df)
}))

combined_sig_count <- rbind(rep_sig_counts, nonrep_sig_counts)

# change order
combined_sig_count$sig <- factor(combined_sig_count$sig, levels = c("rep", "nonrep"))

ggplot(combined_sig_count,
       aes(x = sig, y = sig_count, fill = sig)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.5) +
  theme(axis.text.y   = element_text(size=8),
        axis.text.x   = element_text(size=8),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.line = element_line(colour = "black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("rep" = "indianred3",
                               "nonrep" = "grey90")) +
  ggsignif::geom_signif(comparisons = list(c("rep", "nonrep")),
                        map_signif_level = TRUE,
                        size = 1,
                        margin_top = 0.05,
                        test = "wilcox.test",
                        tip_length = 0.03) +
  labs(x = "RBP shRNA KD experiments in ENCODE",
       y = "Number of SE Events with abs(delPSI) > 0.1") +
  scale_y_continuous(limits = c(0, 2500),
                     breaks = seq(0, 5) * 500,
                     expand = c(0.05, 0))+
  scale_x_discrete(labels = c("rep" = "Rank <= 10 and\nFWER < 0.01",
                              "nonrep" = "Rank > 10 and\nFWER > 0.01")) +
  theme(legend.position = "none")


# manually calculate p-value for manuscript
wilcox.test(subset(combined_sig_count, sig == "rep")$sig_count, subset(combined_sig_count, sig == "nonrep")$sig_count)$p.value


# manually calculate median values
combined_sig_count %>% group_by(sig) %>% summarize(med = mean(sig_count))
