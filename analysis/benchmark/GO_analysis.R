library(enrichR)

num_terms = 5

# get the raw df
raw_benchmark_df <- extract_benchmark_df(sig_weighted_info$raw_benchmark_list[[6]])

rep_count <- table(subset(raw_benchmark_df, gsea_rank <= 10 & gsea_padj < 0.01)$rbp_s2)
non_rep_count <- table(subset(raw_benchmark_df, gsea_rank >= 10 & gsea_padj > 0.01)$rbp_s2)


# rep_count <- table(subset(combined_replicate_result, spl_type == "SE" & sig == "sig")$rbp)  # these are RBPs in top 10 for both HepG2 and K562
# non_rep_count <- table(subset(combined_replicate_result, spl_type == "SE" & sig != "sig")$rbp)  # these are RBPs not in top 10 for both HepG2 and K562

rep_genes <- names(rep_count)[rep_count == 2]
non_rep_genes <- names(non_rep_count)[non_rep_count == 2]

# query GO BP
rep_genes_go <- enrichR::enrichr(rep_genes,
                                 c("GO_Biological_Process_2021"))
non_rep_genes_go <- enrichR::enrichr(non_rep_genes,
                                     c("GO_Biological_Process_2021"))

# extract the data for plotting
rep_genes_top_go <- trim_GO_terms(rep_genes_go$GO_Biological_Process_2021 %>% top_n(wt = -Adjusted.P.value, n = num_terms) %>% arrange(-Adjusted.P.value))
non_rep_genes_top_go <- trim_GO_terms(non_rep_genes_go$GO_Biological_Process_2021 %>% top_n(wt = -Adjusted.P.value, n = num_terms) %>% arrange(-Adjusted.P.value))

rep_genes_top_go$rep <- "top10"
non_rep_genes_top_go$rep <- "not_top10"

rep_genes_top_go$term_short <- factor(rep_genes_top_go$term_short,
                                      levels = rep_genes_top_go$term_short)
non_rep_genes_top_go$term_short <- factor(paste0(non_rep_genes_top_go$term_short, " "),
                                      levels = paste0(non_rep_genes_top_go$term_short, " "))


go_plot_data <- rbind(rep_genes_top_go, non_rep_genes_top_go)

# generate labeller
rep_label <- c("top10" = sprintf("Rank <= 10 & FWER < 0.01\n(N = %s)", length(rep_genes)),
               "not_top10" = sprintf("Rank >= 10 & FWER > 0.01\n(N = %s)", length(non_rep_genes)))
go_plot_data$rep <- factor(go_plot_data$rep,
                           levels = c("top10", "not_top10"))

# go_plot_data$term_short <- factor(go_plot_data$term_short,
#                                   levels = unique(go_plot_data$term_short[order(-go_plot_data$Adjusted.P.value)]))

# generate plot
ggplot(go_plot_data, aes(x = -log10(Adjusted.P.value), y = term_short,
                         fill = rep)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("top10" = "indianred3",
                               "not_top10" = "grey90")) +
  facet_grid(rep ~ ., scales = "free",
             labeller = labeller(rep = rep_label)) +
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
                                    size = 1)
  ) +
  geom_text(data = go_plot_data,
            aes(label = term_short,
                y = term_short),
            hjust = 0,
            x = 1,
            size = 3) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "none") +
  labs(x = "-log10(FDR)",
       y = sprintf("Top %s GO Terms in Biological Process", num_terms))


