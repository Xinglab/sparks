library(dplyr)
library(ggplot2)
library(SPARKS)

##### FUNCTIONS #####

##### ANALYSIS #####
# define the thresholds - these are the ones tested
input_counts <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100)

# define test parameters
testing_method <- "bonferroni"
sig_threshold <- 0.01
rank_threshold <- 10  # looking at top 10 out of 210
sig_weighted_dir <- "/Users/harryyang/research/Xing_Lab/sparks/data/benchmark_final"

# import the benchmark data across the thresholds in 3 splice types
sig_weighted_info_a3ss <- import_benchmark_rerun_data_with_num_events(sig_weighted_dir,
                                                                      input_counts,
                                                                      study = "Final_analysis",
                                                                      testing_method = testing_method,
                                                                      sig_threshold = sig_threshold,
                                                                      rank_threshold = rank_threshold,
                                                                      spl_type = "A3SS")
sig_weighted_info_a5ss <- import_benchmark_rerun_data_with_num_events(sig_weighted_dir,
                                                                      input_counts,
                                                                      study = "Final_analysis",
                                                                      testing_method = testing_method,
                                                                      sig_threshold = sig_threshold,
                                                                      rank_threshold = rank_threshold,
                                                                      spl_type = "A5SS")
sig_weighted_info_se <- import_benchmark_rerun_data_with_num_events(sig_weighted_dir,
                                                                    input_counts,
                                                                    study = "Final_analysis",
                                                                    testing_method = testing_method,
                                                                    sig_threshold = sig_threshold,
                                                                    rank_threshold = rank_threshold,
                                                                    spl_type = "SE")
##### BENCHMARK PLOT OVER THRESHOLDS #####
# generate benchmark plots - w500 h600
generate_benchmark_plot_over_threshold(sig_weighted_info_se)
generate_benchmark_plot_over_threshold(sig_weighted_info_a3ss)
generate_benchmark_plot_over_threshold(sig_weighted_info_a5ss)



##### TOP 10 CDF plot #####
generate_benchmark_cdf_plot_all_sig(extract_benchmark_df(sig_weighted_info_se$raw_benchmark_list[[6]]),
                                    sig_test_method = testing_method,
                                    sig_threshold = sig_threshold)
generate_benchmark_cdf_plot_all_sig(extract_benchmark_df(sig_weighted_info_a5ss$raw_benchmark_list[[6]]),
                                    sig_test_method = testing_method,
                                    sig_threshold = sig_threshold)+ theme(legend.box = "horizontal", legend.position = "bottom", legend.direction = "vertical")
generate_benchmark_cdf_plot_all_sig(extract_benchmark_df(sig_weighted_info_a3ss$raw_benchmark_list[[6]]),
                                    sig_test_method = testing_method,
                                    sig_threshold = sig_threshold)+ theme(legend.box = "horizontal", legend.position = "bottom", legend.direction = "vertical")

#### plot at 30% for all of them ####
top10_se_df <- subset(sig_weighted_info_se$top10_counts,
                      method == "SPARKS" & endsWith(test, '30'))
top10_a3ss_df <- subset(sig_weighted_info_a3ss$top10_counts,
                        method == "SPARKS" & endsWith(test, '30'))
top10_a5ss_df <- subset(sig_weighted_info_a5ss$top10_counts,
                        method == "SPARKS" & endsWith(test, '30'))

top10_se_df$spl_type <- "SE"
top10_a3ss_df$spl_type <- "A3SS"
top10_a5ss_df$spl_type <- "A5SS"

# combine the data for bar plot
top10_plot_df <- do.call(rbind, list(top10_se_df,
                                     top10_a3ss_df,
                                     top10_a5ss_df))

top10_plot_df$spl_type <- factor(top10_plot_df$spl_type,
                                 levels = c("SE", "A5SS", "A3SS"))

ggplot(top10_plot_df) +
  geom_bar(aes(x = spl_type,
               y = value),
           stat = "identity",
           fill = "indianred3",
           color = "black") +
  facet_wrap(~study) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  scale_y_continuous(limits = c(0, 100),
                     expand = c(0, 0)) +
  geom_text(aes(x = spl_type,
                y = value + 3,
                label = value),
            size = 3, color = "indianred3") +
  labs(x = "Splicing Type",
       y = "# Perturbed RBP with\nRank <= 10 & FWER < 0.01") # 400 x 300


##### P-VALUE CURVE OVER COUNT CHANGE #####
# import p-value data for different counts
es_pval_plot_df <- do.call(rbind, lapply(seq(length(input_counts)), function(idx){
  print(input_counts[idx])

  # query the benchmark df
  test_benchmark_df <- sig_weighted_info$raw_benchmark_list[[idx]]$extra

  # mutate the data as necessary
  es_diff_df <- test_benchmark_df[, c("S1", "S2", "gsea_combined_pval")]

  es_diff_plot_df <- reshape2::melt(es_diff_df)
  es_diff_plot_df$count <- input_counts[idx]
  return(es_diff_plot_df)
}))

# factorize for manual coloring
es_pval_plot_df$count <- factor(es_pval_plot_df$count, levels = input_counts)

# assign colors
color_list <- colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(15)[4:15]
names(color_list) <- input_counts

# generate score diff plot
ggplot(es_pval_plot_df) +
  stat_ecdf(aes(x = -log10(value), group = count, color = count), alpha = 1, size = 1) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = color_list) +
  labs(x = "-log10(p-value)",
       y = "CDF for each RBP KD dataset",
       color = "Filter Threshold") +
  coord_cartesian(xlim = c(0, 20))


##### CDF PLOT FOR OVERLAP #####
kd_library_all <- readRDS("/Users/harryyang/transfer/Clean_minfiltered.A5SS.library.rds")

# distribution test
# query null events
event_list <- list()
# extract all sig events in the library
dummy <- lapply(names(kd_library_all), function(experiment){
  sig_events <- unlist(extract_GSEA_significant_events(kd_library_all[[experiment]]))
  event_list[[length(event_list) + 1]] <<- sig_events
  return()
})

# geenrate list and weight for them
null_event <- unique(unlist(event_list))
null_weight <- table(unlist(event_list))

null_weight_plot_df <- as.data.frame(null_weight)
colnames(null_weight_plot_df) <- c("event", "freq")

ggplot(null_weight_plot_df, aes(x = freq)) +
  stat_ecdf() +
  scale_y_continuous(breaks = seq(0, 4)/4,
                     labels = round(seq(0, 4) / 4 * length(null_event), 0),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
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
  labs(x = "# RBP KD experiments in the library\nthat an AS event is significant",
       y = sprintf("Cumulative # AS events significant in the library\n(Total N = %s)", length(null_event))) +
  scale_x_continuous(limits = c(0, max(null_weight_plot_df$freq)),
                     expand = c(0, 0),
                     breaks = round(c(0, median(null_weight_plot_df$freq), 100, 200, 300, max(null_weight_plot_df$freq)))) +
  # add median
  geom_vline(xintercept = median(null_weight_plot_df$freq),
             linetype = "dashed",
             color = "indianred1") +
  geom_hline(yintercept = 0.5 ,
             linetype = "dashed",
             color = "indianred1")



##### SAME RBP RANK DISTRIBUTION #####
benchmark_rank_score_combined <- do.call(rbind, lapply(seq(11), function(idx){

  # get the rank distribution
  benchmark_rank_scores <- data.frame(rank = subset(sig_weighted_info$raw_benchmark_list[[idx]]$extra, rbp_s1 == rbp_s2)$gsea_rank,
                                      count = input_counts[idx])

  print(head(benchmark_rank_scores))
  return(benchmark_rank_scores)
}))

benchmark_rank_score_combined$count <- factor(benchmark_rank_score_combined$count,
                                              levels = input_counts)
# assign colors
color_list <- colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(15)[4:15]
names(color_list) <- input_counts

ggplot(benchmark_rank_score_combined, aes(x = rank,
                                          color = count,
                                          group = count)) +
  geom_density()+
  geom_vline(xintercept = 10,
             linetype = "dashed",
             color = "indianred2", size = 0.5) +
  labs(x = "SPARKS Enrichment Score\nfor RBP pairs with same sign between ES+ and ES-",
       y = "Density") +
  scale_color_brewer(palette = "Paired") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))


##### BENCHMARK FIGURE #####
benchmark_count_df <- sig_weighted_info_a3ss$top10_counts

benchmark_count_df$param <- unlist(lapply(benchmark_count_df$test,
                                          function(x) strsplit(as.character(x), "= ")[[1]][2]))
# set manual order to treat it like a character
benchmark_count_df$param <- factor(benchmark_count_df$param,
                                   levels = unique(benchmark_count_df$param))

# map it back to the original study
var_names <- c("Pearson's Rho-based Method", "Pearson's Rho-based Method",
               "Kendall's Tau-based Method", "Kendall's Tau-based Method",
               "SPARKS", "SPARKS")
names(var_names) <- unique(benchmark_count_df$variable)


benchmark_count_df$original_method <- var_names[benchmark_count_df$variable]

# seperate sig vs. ns
benchmark_count_df$sig <- ifelse(endsWith(as.character(benchmark_count_df$variable), 'ns'), 'Notsig', 'Sig')

benchmark_count_plot_df <- benchmark_count_df %>% group_by(method, param, original_method, sig) %>% summarise(value = sum(value))
benchmark_count_plot_df$original_method <- factor(benchmark_count_plot_df$original_method,
                                                  levels = c("SPARKS",
                                                             "Pearson's Rho-based Method",
                                                             "Kendall's Tau-based Method"))
benchmark_count_plot_df$sig <- factor(benchmark_count_plot_df$sig,
                                      levels = c('Sig', 'Notsig'))

ggplot(benchmark_count_plot_df,
       aes(x = param, y = value,
           color = original_method,
           # group = sig,
           group = interaction(original_method, sig),
           shape = sig)) +
  geom_point() + geom_line(aes(linetype = sig))+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  scale_linetype_discrete(labels = c("Sig" = "FWER < 0.01",
                                     "Notsig" = "FWER > 0.01")) +
  scale_shape_discrete(labels = c("Sig" = "FWER < 0.01",
                                  "Notsig" = "FWER > 0.01")) +
  labs(x = "Overlap Percentage Threshold",
       y = "# RBP experiments with\nthe perturbed RBP ranked within 10\nacross HepG2 and K562 (N = 420)",
       color = "Method",
       linetype = "Significance",
       shape = "Significance")


##### SRSF1 example figure #####
generate_strip_lollipop_plot_vertical(list(SE = subset(sig_weighted_info_se$raw_benchmark_list[[6]]$extra, S2 == "K562_SRSF1_shRNA")),
                                      spl_type = 'SE',
                                      rbp_of_interest = c("SRSF1"),
                                      sig_test = T,
                                      sig_test_method = "bonferroni",
                                      max_score_threshold = 1.5) # 350 x 200
generate_strip_lollipop_plot_vertical(list(SE = subset(sig_weighted_info_se$raw_benchmark_list[[6]]$extra, S2 == "HepG2_SRSF1_shRNA")),
                                      spl_type = 'SE',
                                      rbp_of_interest = c("SRSF1"),
                                      sig_test = T,
                                      sig_test_method = "bonferroni",
                                      max_score_threshold = 1.5)


##### HEATMAP FOR CORR #####
raw_benchmark_data <- sig_weighted_info$raw_benchmark_list[[6]]$extra

test_cell <- 'HepG2'
other_cell <- 'K562'
test_cell <- "K562"
other_cell <- "HepG2"

# take subset for plotting
plot_bench_df <- subset(raw_benchmark_data, startsWith(S1, other_cell) & startsWith(S2, test_cell))

test_rbp_list <- data.table::fread('/Users/harryyang/research/Xing_Lab/START/other_data/SF_list.txt', header = F)$V1

plot_bench_df <- subset(plot_bench_df, !(rbp_s1 %in% test_rbp_list) & !(rbp_s2 %in% test_rbp_list))
# plot_bench_df <- subset(plot_bench_df, rbp_s1 %in% test_rbp_list & rbp_s2 %in% test_rbp_list)

# annotate frequent ones
most_sig_rbps <- (plot_bench_df %>% group_by(rbp_s1) %>% summarize(total_sig = sum(gsea_rank <= 10)) %>% filter(total_sig > 21))$rbp_s1

most_sig_label <- ifelse(plot_bench_df$rbp_s1 %in% most_sig_rbps, plot_bench_df$rbp_s1, "")

ggplot(plot_bench_df,
       aes(x = rbp_s2,
           y = rbp_s1,
           fill = gsea_rank)) +
  geom_raster() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  # theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  # axis.ticks = element_blank()) +
  # scale_y_discrete)) +
  scale_fill_distiller(palette = "RdPu",
                       direction = -1,
                       limits = c(1, 10),
                       na.value = "white",
                       breaks = c(1, 5, 10)) +
  labs(x = sprintf("RBP KD AS Profiles from %s",
                   test_cell),
       y = sprintf("RBP KD AS Signatures from %s",
                   other_cell),
       fill = "SPARKS Rank") +
  geom_abline(intercept = c(0, 0), slope = 1, size = 0.25, linetype = "dashed") +
  theme(axis.text = element_blank())

# score plot
ggplot(plot_bench_df,
       aes(x = rbp_s2,
           y = rbp_s1,
           fill = gsea_score)) +
  geom_raster() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  # scale_y_discrete(label = most_sig_label) +
  scale_fill_distiller(palette = "RdBu",
                       direction = -1,
                       limits = c(-2, 2)) +
  labs(x = sprintf("RBP KD AS Profiles from %s",
                   test_cell),
       y = sprintf("RBP KD AS Signatures from %s",
                   other_cell),
       fill = "SPARKS ES")




##### RANK DISTRIBUTION TEST #####
raw_benchmark_df <- sig_weighted_info_se$raw_benchmark_list[[6]]$extra

# update the rbp s1 for same vs. not
raw_benchmark_df$new_rbp <- ifelse(raw_benchmark_df$rbp_s1 == raw_benchmark_df$rbp_s2,
                                   "same",
                                   raw_benchmark_df$rbp_s1)

# update order
raw_benchmark_df$new_rbp <- factor(raw_benchmark_df$new_rbp,
                                   levels = c(unique(raw_benchmark_df$rbp_s1), "same"))

# calculate p-value for the difference
perturbed_rbp_dist <- subset(raw_benchmark_df, new_rbp == "same")$gsea_rank

pval_result <- do.call(rbind, lapply(unique(raw_benchmark_df$rbp_s1), function(test_rbp){
  other_rbp_dist <- subset(raw_benchmark_df, new_rbp == test_rbp)$gsea_rank

  ks_pval <- ks.test(perturbed_rbp_dist, other_rbp_dist)$p.value
  wilcox_pval <- wilcox.test(perturbed_rbp_dist, other_rbp_dist)$p.value

  pval_result_df <- data.frame(rbp = test_rbp,
                               ks = ks_pval,
                               wilcox_pval = wilcox_pval)
  return(pval_result_df)
}))

median_p <- median(pval_result$wilcox_pval)
median_p_data <- data.frame(p = median_p,
                            x = 180,
                            y = 0.013)

ggplot(raw_benchmark_df) +
  geom_density(aes(x = gsea_rank,
                   color = ifelse(new_rbp == 'same', 'same', 'other'),
                   group = new_rbp),
               bounds = c(1, 210)) +
  scale_color_manual(values = c('same' = 'indianred2',
                                'other' = 'grey85'),
                     label = c('same' = "Pertrubed RBP",
                               'other' = "Other RBPs (n = 210)")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) +
  theme(legend.position = c(0.55, 0.75),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  labs(x = "SPARKS Rank",
       y = "Density",
       color = "RBP") +
  scale_y_continuous(expand = c(0.05, 0)) +
  geom_text(data = median_p_data,
            aes(label = sprintf("median p < %s",
                                signif(p, 2))),
            size = 3, x = 80, y = 0.024,
            color = "indianred2",
            fontface = "bold",
            hjust = 0.5) # w 300 x h 400

# calculate the stats
raw_benchmark_df %>% filter(new_rbp == "same") %>% summarize(mm = median(gsea_rank))
mean((raw_benchmark_df %>% filter(new_rbp != "same") %>% group_by(new_rbp) %>% summarize(mm = median(gsea_rank)))$mm)

### PVALUE
pval_plot_data <- raw_benchmark_df %>% group_by(new_rbp) %>% summarize(total_sig = sum(gsea_padj < 0.01))
pval_plot_data$filtered_rbp <- ifelse(pval_plot_data$total_sig > 125,
                                      as.character(pval_plot_data$new_rbp), NA)
pval_plot_data$filtered_rbp[pval_plot_data$filtered_rbp == "same"] <- "Perturbed RBP"

ggplot(pval_plot_data,
       aes(x = "",
           fill = ifelse(new_rbp == 'same', 'same', 'other'),
           # group = ifelse(new_rbp == 'same', 'same', 'other'),
           y = total_sig)) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5, binwidth = 5) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c('same' = 'indianred2',
                               'other' = 'grey75'),
                    label = c('same' = "Pertrubed RBP",
                              'other' = "Other RBPs (n = 210)")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) +
  labs(x = "RBP perturbed in both\nHepG2 and K562 (n = 210)",
       y = "# Significant experiment\n(out of 420 RBP KD experiments)",
       fill = "RBP") +
  scale_y_continuous(limits = c(0, 200)) +
  theme(legend.position = 'none') +
  ggrepel::geom_label_repel(aes(label = filtered_rbp),
                            size = 2, force =100, direction = "x")


### RANK + PVAL

both_plot_data <- raw_benchmark_df %>% group_by(new_rbp) %>% summarize(total_sig = sum(gsea_padj < 0.01 & gsea_rank <= 10))
both_plot_data$filtered_rbp <- ifelse(both_plot_data$total_sig > 45,
                                      as.character(both_plot_data$new_rbp), NA)
both_plot_data$filtered_rbp[both_plot_data$filtered_rbp == "same"] <- "Perturbed RBP"

ggplot(both_plot_data,
       aes(x = "",
           fill = ifelse(new_rbp == 'same', 'same', 'other'),
           # group = ifelse(new_rbp == 'same', 'same', 'other'),
           y = total_sig)) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1, binwidth = 2) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c('same' = 'indianred2',
                               'other' = 'grey75'),
                    label = c('same' = "Pertrubed RBP",
                              'other' = "Other RBPs (n = 210)")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) +
  labs(x = "RBP perturbed in both\nHepG2 and K562 (n = 210)",
       y = "# Experiments with Rank <= 10 & FWER < 0.01\n(out of 420 RBP KD experiments)",
       fill = "RBP") +
  scale_y_continuous(limits = c(0, 160),
                     breaks = seq(0, 4) * 40) +
  theme(legend.position = 'none') +
  ggrepel::geom_label_repel(aes(label = filtered_rbp),
                            size = 2.5, force = 30, direction = "both",
                            # ylim = c(50, 150),
                            nudge_x = 0.5)




test_rbp_list <- unique(raw_benchmark_df$rbp_s1)

pval_result_cross <- do.call(rbind, lapply(test_rbp_list, function(test_rbp){
  print(test_rbp)
  pval_result_cross2 <- do.call(rbind, lapply(test_rbp_list[test_rbp_list != test_rbp], function(test_rbp2){
    print(test_rbp2)
    other_rbp_dist <- subset(raw_benchmark_df, new_rbp == test_rbp)$gsea_rank
    other_rbp_dist2 <- subset(raw_benchmark_df, new_rbp == test_rbp2)$gsea_rank

    ks_pval <- ks.test(other_rbp_dist2, other_rbp_dist)$p.value
    wilcox_pval <- wilcox.test(other_rbp_dist2, other_rbp_dist)$p.value

    pval_result_df <- data.frame(rbp1 = test_rbp,
                                 rbp2 = test_rbp2,
                                 ks = ks_pval,
                                 wilcox_pval = wilcox_pval)
    return(pval_result_df)
  }))
  return(pval_result_cross2)
}))

# generate plot
same_rbp_test_df <- data.frame(p = pval_result$wilcox_pval,
                               test = "same")
other_rbp_test_df <- data.frame(p = pval_result_cross$wilcox_pval,
                                test = "other")
combined_plot_df <- rbind(same_rbp_test_df, other_rbp_test_df)

ggplot(combined_plot_df, aes(x = test, y = -log10(p))) +
  geom_jitter(alpha = 0.1, width = 0.25) +
  geom_boxplot(width = 0.5, alpha = 0.5) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) +
  labs(x = "Status",
       y = "-log10(p-value)") +
  scale_x_discrete(labels = c("same" = "Perturbed RBP vs. Other RBP",
                              "other" = "Other RBP vs. Other RBP"))





##### SAME VS OTHER CELL LINE TEST ####
same_benchmark_df <- sig_weighted_info_se$raw_benchmark_list[[6]]$intra
other_benchmark_df <- sig_weighted_info_se$raw_benchmark_list[[6]]$extra

# query score
study_cell_line <- "K562"
study_rbp <- "SRSF1"
test_rbp <- "PTBP1"


study_rbps <- unique(other_benchmark_df$rbp_s1)

combined_score_df <- do.call(rbind, lapply(study_rbps, function(study_rbp){
  # gather significant rbp experiments
  sig_benchmark_df <- subset(other_benchmark_df, S2 == sprintf("%s_%s_shRNA", study_cell_line, study_rbp) & gsea_padj < 0.0000001)

  test_rbps <- unique(sig_benchmark_df$rbp_s1)


  score_sum_df <- do.call(rbind, lapply(test_rbps[!(test_rbps == study_rbp)], function(test_rbp){
    print(test_rbp)

    # query the score
    same_df <- subset(same_benchmark_df, S2 == sprintf("%s_%s_shRNA", study_cell_line, study_rbp) &
                        rbp_s1 ==  test_rbp)
    other_df <- subset(other_benchmark_df, S2 == sprintf("%s_%s_shRNA", study_cell_line, study_rbp) &
                         rbp_s1 ==  test_rbp)
    same_score <- abs(same_df$gsea_pos_score) + abs(same_df$gsea_neg_score)
    other_score <- abs(other_df$gsea_pos_score) + abs(other_df$gsea_neg_score)

    score_df <- data.frame(
      study_rbp = study_rbp,
      test_rbp = test_rbp,
      same_score = same_score,
      other_score = other_score
    )
    return(score_df)
  }))
  print(score_sum_df)
  # score_sum_df$diff <- abs(score_sum_df$same_score) - abs(score_sum_df$other_score)

  return(score_sum_df)
}))

combined_score_df$diff <- combined_score_df$same_score - combined_score_df$other_score

# calculate Num sig
combined_score_df <- do.call(rbind, lapply(study_rbps, function(study_rbp){
  cl_score_df <- do.call(rbind, lapply(c("K562", "HepG2"), function(study_cell_line){


    print(study_rbp)

    # query the score
    same_score <- dim(subset(same_benchmark_df, S2 == sprintf("%s_%s_shRNA", study_cell_line, study_rbp) & gsea_padj < 0.01))[1]
    other_score <- dim(subset(other_benchmark_df, S2 == sprintf("%s_%s_shRNA", study_cell_line, study_rbp) & gsea_padj < 0.01))[1]

    score_df <- data.frame(
      study_rbp = study_rbp,
      study_cell_line = study_cell_line,
      same_score = same_score,
      other_score = other_score
    )

    return(score_df)}))
  return(cl_score_df)
}))

# combined_score_df$diff <- combined_score_df$same_score - combined_score_df$other_score

combined_score_melt <- reshape2::melt(combined_score_df)
combined_score_melt$variable <- factor(combined_score_melt$variable,
                                       levels = c("other_score", "same_score"))

ggplot(combined_score_melt, aes(x = variable, y = value, group = study_rbp)) +
  geom_point(alpha = 0.3) +
  geom_line(alpha = 0.3) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) +
  labs(x = "Testing Condition",
       y = "# Siginificant experiment (out of 210)") +
  scale_x_discrete(labels = c("same_score" = "Same cell line",
                              "other_score" = "Other cell line")) +
  scale_y_continuous(limits = c(0, 200)) +
  facet_wrap(~study_cell_line)

