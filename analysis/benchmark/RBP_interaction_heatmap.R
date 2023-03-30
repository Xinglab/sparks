##### SAME VS OTHER CELL LINE TEST ####
same_benchmark_df <- sig_weighted_info_se$raw_benchmark_list[[6]]$intra
other_benchmark_df <- sig_weighted_info_se$raw_benchmark_list[[6]]$extra

# query list ofRBPs
study_rbps <- unique(other_benchmark_df$rbp_s1)

# subset significant experiments
same_sig_subset <- subset(same_benchmark_df, gsea_rank <= 10 & gsea_padj < 0.01)
other_sig_subset <- subset(other_benchmark_df, gsea_rank <= 10 & gsea_padj < 0.01)

# count occurrences separately - both for other, one for same?
other_sig_counts <- subset(as.data.frame(table(other_sig_subset[, c("rbp_s1", "rbp_s2")])), Freq >= 2)
same_sig_counts <- subset(as.data.frame(table(same_sig_subset[, c("rbp_s1", "rbp_s2")])), Freq >= 2)

#

sig_counts <- intersect(other_sig_counts[, c("rbp_s1", "rbp_s2")], same_sig_counts[, c("rbp_s1", "rbp_s2")])

sig_rbps <- union(sig_counts$rbp_s1, sig_counts$rbp_s2)

# average the score
raw_score_df <- subset(other_benchmark_df, rbp_s1 %in% sig_rbps & rbp_s2 %in% sig_rbps)
raw_score_df <- rbind(subset(other_benchmark_df, rbp_s1 %in% sig_rbps & rbp_s2 %in% sig_rbps)[, c("rbp_s1", "rbp_s2", "gsea_score")],
                      subset(same_benchmark_df, rbp_s1 %in% sig_rbps & rbp_s2 %in% sig_rbps)[, c("rbp_s1", "rbp_s2", "gsea_score")])

# add the reverse for symmetric average
raw_score_df_copy <- raw_score_df
raw_score_df_copy$rbp_s1 <- raw_score_df$rbp_s2
raw_score_df_copy$rbp_s2 <- raw_score_df$rbp_s1
combined_score_df <- rbind(raw_score_df, raw_score_df_copy)

# mutate the data as necessary
plot_bench_df <-  combined_score_df %>% group_by(rbp_s1, rbp_s2) %>% summarize(mean_score = mean(gsea_score))
plot_bench_df$rbp_s1 <- factor(plot_bench_df$rbp_s1,
                               levels = rev(sort(unique(plot_bench_df$rbp_s1))))

# add points to the sig interactions
sig_ixns <- sig_counts
sig_ixns$Freq <- "*"
sig_ixns_cp <- sig_ixns
sig_ixns_cp$rbp_s1 <- sig_ixns$rbp_s2
sig_ixns_cp$rbp_s2 <- sig_ixns$rbp_s1

# combined_sig_ixns <- rbind(sig_ixns, sig_ixns_cp)
combined_sig_ixns_bidirectional <- combined_sig_ixns %>% group_by(rbp_s1, rbp_s2) %>% summarize(a = n())  # test code for bidirectional vs. single directional test



ggplot(plot_bench_df) +
  geom_raster(aes(x = rbp_s2,
                  y = rbp_s1,
                  fill = mean_score)) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  scale_fill_distiller(palette = "RdBu",
                       direction = -1,
                       limits = c(-max(plot_bench_df$mean_score), max(plot_bench_df$mean_score))) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  geom_text(data = combined_sig_ixns,
            aes( x = rbp_s2, y = rbp_s1, label = Freq), vjust = 0.75) +
  labs(x = "RBPs",
       y = "RBPs",
       fill = "Mean ES")



# clustering for ordering
aaa <- reshape2::dcast(plot_bench_df, rbp_s1 ~ rbp_s2, value.var = 'mean_score')
bbb <- hclust(dist(aaa, method = "manhattan"))
rbp_order <- aaa$rbp_s1[bbb$order]

# ordering based on clustering
plot_bench_df_cp <- plot_bench_df
plot_bench_df_cp$rbp_s1 <- factor(plot_bench_df_cp$rbp_s1,
                                  levels = rbp_order)
plot_bench_df_cp$rbp_s2 <- factor(plot_bench_df_cp$rbp_s2,
                                  levels = rbp_order)

# generate plot
ggplot(plot_bench_df_cp) +
  geom_raster(aes(x = rbp_s2,
                  y = rbp_s1,
                  fill = mean_score)) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  scale_fill_distiller(palette = "RdBu",
                       direction = -1,
                       limits = c(-max(plot_bench_df$mean_score), max(plot_bench_df$mean_score))) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  geom_text(data = combined_sig_ixns,
            aes( x = rbp_s2, y = rbp_s1, label = Freq), vjust = 0.75) +
  labs(x = "RBPs",
       y = "RBPs",
       fill = "Mean ES")

