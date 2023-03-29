test_rbp_list <- data.table::fread('/Users/harryyang/research/Xing_Lab/START/other_data/SF_list.txt', header = F)$V1

bench_plot_df <- subset(sig_weighted_info_se$raw_benchmark_list[[6]]$extra, rbp_s1 == rbp_s2)

test_sig_threshold <- 0.01
test_rank <- 10

bench_plot_df$status <- ifelse(bench_plot_df$gsea_rank <= test_rank & bench_plot_df$gsea_padj <= test_sig_threshold,
                               "sig", "notsig")
# reorder
bench_plot_df$status <- factor(bench_plot_df$status,
                               levels = c("notsig", "sig"))



bench_plot_df$sf_status <- ifelse(bench_plot_df$rbp_s2 %in% test_rbp_list,
                                  "SF", "Other RBP")
# reorder
bench_plot_df$sf_status <- factor(bench_plot_df$sf_status,
                                  levels = c("SF", "Other RBP"))

# change name for x axis
num_sf <- length(intersect(unique(bench_plot_df$rbp_s2), test_rbp_list))
sf_status_names <- c(sprintf("Splicing Factor\n(n = %s)", num_sf),
                     sprintf("Other RBPs\n(n = %s)", 210 - num_sf))
names(sf_status_names) <- c("SF", "Other RBP")

ggplot(bench_plot_df,
       aes(x = sf_status)) +
  geom_bar(aes(fill = status), stat = "count", position = "fill") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) +
  scale_fill_manual(values = c("notsig" = "lavenderblush3",
                               "sig" = "indianred3"),
                    labels = c("notsig" = "Rank > 10 or FWER > 0.01",
                               "sig" = "Rank <= 10 and FWER <= 0.01"))+
  scale_x_discrete(label = sf_status_names) +
  labs(x = "Known Splicing Factors",
       y = "Proportion",
       fill = "SPARKS Result on the Perturbed RBP") +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  labs(subtitle = sprintf("p = %s", signif(fisher.test(table(bench_plot_df[, c("status", "sf_status")]))$p.value, 3))) +
  theme(plot.subtitle = element_text(size = 10))

# proportion for manuscript
(bench_plot_df) %>% group_by(sf_status) %>% summarize(mean = sum(status == 'sig') / n())

# p-value for manuscript
fisher.test(table(bench_plot_df[, c("status", "sf_status")]))$p.value



### threshold testing

test_rank_thresholds <- c(10, 20, 30, 40, 50, 75, 105)
test_sig_thresholds <- c(0.01, 0.05, 0.1, 0.2)

test_df <- bench_plot_df

#
threshold_test_df <- do.call(rbind, lapply(test_rank_thresholds, function(test_rank_threshold){
  sig_test_df <- do.call(rbind, lapply(test_sig_thresholds, function(test_sig_threshold){
    test_df$status <- ifelse(bench_plot_df$gsea_rank <= test_rank_threshold & bench_plot_df$gsea_padj <= test_sig_threshold,
                             "sig", "notsig")
    test_result_df <- as.data.frame(table(test_df[, c("status", "sf_status")]))
    test_result_df$rank_threshold <- test_rank_threshold
    test_result_df$sig_threshold <- test_sig_threshold
    print(test_result_df)
    return(test_result_df)
  }))
  return(sig_test_df)
}))

rank_plot_df <- subset(threshold_test_df, sig_threshold == 0.01)


ggplot(rank_plot_df,
       aes(x = sf_status)) +
  geom_bar(aes(fill = status, y = Freq), stat = "identity", position = "fill") +
  facet_wrap(~ rank_threshold)+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) +
  scale_fill_manual(values = c("notsig" = "lavenderblush3",
                               "sig" = "indianred3"),
                    labels = c("notsig" = "Rank > 10 or FWER > 0.01",
                               "sig" = "Rank <= 10 and FWER <= 0.01"))+
  scale_x_discrete(label = sf_status_names) +
  labs(x = "Known Splicing Factors",
       y = "Proportion",
       fill = "SPARKS Result on the Perturbed RBP") +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))


padj_plot_df <- subset(threshold_test_df, rank_threshold == 10)


ggplot(threshold_test_df,
       aes(x = sf_status)) +
  geom_bar(aes(fill = status, y = Freq), stat = "identity", position = "fill") +
  facet_grid(rank_threshold ~ sig_threshold)+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) +
  scale_fill_manual(values = c("notsig" = "lavenderblush3",
                               "sig" = "indianred3"),
                    labels = c("notsig" = "Rank > Rank Threshold or \nFWER > Significance Threshold",
                               "sig" = "Rank <= Rank Threshold and \nFWER <= Significance Threshold"))+
  scale_x_discrete(label = sf_status_names) +
  labs(x = "Known Splicing Factors",
       y = "Proportion",
       fill = "SPARKS Result on the Perturbed RBP") +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  geom_text(data = subset(threshold_test_df %>% group_by(rank_threshold, sig_threshold, sf_status) %>% mutate(prop = Freq / sum(Freq)),
                          status == 'sig'),
            aes(x = sf_status, y = prop + 0.15, label = Freq),
            size = 3, color = "indianred3") +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  # 400 x 600


##### OLD CODE #####
# define the criteria for color
bench_plot_df$test_new <- sprintf("Rank %s 10 & %s %s %s",
                                  ifelse(bench_plot_df$gsea_rank <= 10,
                                         "<=", ">"),
                                  ifelse(testing_method == "BH",
                                         "FDR", "FWER"),
                                  ifelse(bench_plot_df$gsea_padj <= sig_threshold,
                                         "<", ">="),
                                  sig_threshold)
# reorder
bench_plot_df$test_new <- factor(bench_plot_df$test_new,
                                 levels = c("Rank > 10 & FWER >= 0.01",
                                            "Rank > 10 & FWER < 0.01",
                                            "Rank <= 10 & FWER >= 0.01",
                                            "Rank <= 10 & FWER < 0.01"))



bench_plot_df$sf_status <- ifelse(bench_plot_df$rbp_s2 %in% test_rbp_list,
                                  "SF", "Other RBP")
# reorder
bench_plot_df$sf_status <- factor(bench_plot_df$sf_status,
                                  levels = c("SF", "Other RBP"))

# change name for x axis
num_sf <- length(intersect(unique(bench_plot_df$rbp_s2), test_rbp_list))
sf_status_names <- c(sprintf("Splicing Factor\n(n = %s)", num_sf),
                     sprintf("Other RBPs\n(n = %s)", 210 - num_sf))
names(sf_status_names) <- c("SF", "Other RBP")

ggplot(bench_plot_df,
       aes(x = sf_status)) +
  geom_bar(aes(fill = test_new), stat = "count", position = "fill") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) +
  scale_fill_manual(values = c("Rank > 10 & FWER >= 0.01" = "lavenderblush3",
                               "Rank > 10 & FWER < 0.01" = "skyblue",
                               "Rank <= 10 & FWER >= 0.01" = "lightcoral",
                               "Rank <= 10 & FWER < 0.01" = "orchid3"))+
  scale_x_discrete(label = sf_status_names) +
  labs(x = "Known Splicing Factors",
       y = "Proportion",
       fill = "SPARKS Result on the Perturbed RBP") +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))



