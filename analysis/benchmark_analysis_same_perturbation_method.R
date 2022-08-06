##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(pbmcapply)

##### FUNCTIONS #####
# calculate rank
calculate_benchmark_rank <- function(result_df){

  # process each target study result separately
  target_study_list <- unique(result_df$S1)

  rank_anno_result_df <- do.call(rbind, lapply(target_study_list, function(target_study){
    target_result_df <- subset(result_df, S1 == target_study)


    # calculate rank
    target_result_df$lincor_rank <- rank(-target_result_df$score)
    target_result_df$concord_rank <- rank(-target_result_df$concordance)
    target_result_df$gsea_rank <- rank(-target_result_df$gsea_score)
    return(target_result_df)
  }))

  return(rank_anno_result_df)
}


# extract results within replicate
extract_benchmark_replicate_result <- function(result_df){
  # determine rank for the same ones
  rep_result <- subset(result_df, rbp_s1 == rbp_s2)# & abs(num_eff) > 10)

  # generate plot
  rep_plot_df <- rep_result[, c("lincor_rank", "concord_rank", "gsea_rank", "rbp_s1")]
  return(rep_plot_df)
}


# calculate CDF
calculate_CDF_rank_for_benchmark <- function(result_df){
  rep_plot_df <- extract_benchmark_replicate_result(result_df)

  cdf_df <- do.call(rbind, lapply(c("lincor_rank", "concord_rank", "gsea_rank"), function(x) {

    cumsum_df <- as.data.frame(unlist(lapply(seq(unique(result_df$S1)), function(y)
      return(sum(rep_plot_df[, x] <= y ))
    )))

    colnames(cumsum_df) <- c("cum_count")
    res_df <- data.frame(rank = as.numeric(rownames(cumsum_df)),
                         count = cumsum_df$cum_count,
                         method = x)

    return(res_df)
  }))
  return(cdf_df)
}



##### ANALYSIS #####
## process cluster-run benchmark data
# read in benchmark result data

sf_list <- data.table::fread('/Users/harryyang/research/Xing_Lab/START/other_data/SF_list.txt', header = F)$V1


count_threshold <- 20
test_method <- "average"

spl_types <- c("SE", "A5SS", "A3SS")

hepg2_k562_benchmark_result_list <- list()
k562_hepg2_benchmark_result_list <- list()


all_cdf_plot_df <- do.call(rbind, lapply(spl_types, function(spl_type){
  print(spl_type)
  # input_file <- sprintf('/Users/harryyang/research/Xing_Lab/sparks/data/benchmark_result.shRNA.%s.df.txt', spl_type)
  input_file <- sprintf('/Users/harryyang/research/Xing_Lab/sparks/data/benchmark_result.CRISPR.%s.df.txt', spl_type)

  benchmark_result_df <- as.data.frame(data.table::fread(input_file))

  # annotate cell lin and RBP for easier processing
  benchmark_result_df$study_s1 <- unlist(lapply(benchmark_result_df$S1, function(x) strsplit(x, "_")[[1]][1]))
  benchmark_result_df$study_s2 <- unlist(lapply(benchmark_result_df$S2, function(x) strsplit(x, "_")[[1]][1]))
  benchmark_result_df$rbp_s1 <- unlist(lapply(benchmark_result_df$S1, function(x) strsplit(x, "_")[[1]][2]))
  benchmark_result_df$rbp_s2 <- unlist(lapply(benchmark_result_df$S2, function(x) strsplit(x, "_")[[1]][2]))

  ## calculate rank for each condition
  # S1 = the target study
  # S2 = the library
  hepg2_k562_result_df <- subset(benchmark_result_df, study_s1 == "HepG2" & study_s2 == "K562")
  k562_hepg2_result_df <- subset(benchmark_result_df, study_s1 == "K562" & study_s2 == "HepG2")
  hepg2_hepg2_result_df <- subset(benchmark_result_df, study_s1 == "HepG2" & study_s2 == "HepG2")
  k562_k562_result_df <- subset(benchmark_result_df, study_s1 == "K562" & study_s2 == "K562")

  # calculate benchmark rank for cross cell line
  hepg2_k562_result_df_anno <- calculate_benchmark_rank(hepg2_k562_result_df)
  k562_hepg2_result_df_anno <- calculate_benchmark_rank(k562_hepg2_result_df)

  ## extract the rank for replicate for downstream analysis
  hepg2_k562_rank_df <- extract_benchmark_replicate_result(hepg2_k562_result_df_anno)
  k562_hepg2_rank_df <- extract_benchmark_replicate_result(k562_hepg2_result_df_anno)

  # annotate spl type for data structure
  hepg2_k562_rank_df$spl_type <- spl_type
  k562_hepg2_rank_df$spl_type <- spl_type

  # add it to outside list for analysis
  hepg2_k562_benchmark_result_list[[spl_type]] <<- hepg2_k562_rank_df
  k562_hepg2_benchmark_result_list[[spl_type]] <<- k562_hepg2_rank_df


  ## Generate CDF plot
  # gather data
  hepg2_cdf_plot_df <- calculate_CDF_rank_for_benchmark(hepg2_k562_result_df_anno)
  k562_cdf_plot_df <- calculate_CDF_rank_for_benchmark(k562_hepg2_result_df_anno)

  hepg2_cdf_plot_df$study <- "HepG2 using K562 Library"
  k562_cdf_plot_df$study <- "K562 using HepG2 Library"

  # merge the data
  top_cdf_plot_df <- rbind(subset(hepg2_cdf_plot_df, rank <= 10),
                           subset(k562_cdf_plot_df, rank <= 10))

  top_cdf_plot_df$spl_type <- spl_type
  return(top_cdf_plot_df)
}))

# make label
test_names <- c("lincor_rank", "concord_rank", "gsea_rank")
new_test_names <- c("Pearson's Rho-based Method", "Kendall's Tau-based Method", "SPARKS")
names(new_test_names) <- test_names

# re order the splicing types
all_cdf_plot_df$spl_type <- factor(all_cdf_plot_df$spl_type,
                                   levels = spl_types)

# generate plot
ggplot(all_cdf_plot_df, aes(x = rank, y = count, color = method)) + geom_point() + geom_line() +
  scale_color_discrete(labels = new_test_names) +
  scale_x_continuous(breaks = seq(-1, 10),
                     limits = c(1, 10)) +
  # scale_y_continuous(limits = c(0, 90), breaks = c(seq(0, 6) * 15)) +  # scale for shRNA KD
  scale_y_continuous(limits = c(0, 35), breaks = c(seq(0, 7) * 5)) +  # scale for CRISPR KD
  labs(x = "Rank of the target RBP from each method",
       y = "Cumulative Sum of RBP KD Experiments",
       color = "Method") +
  facet_grid(study ~ spl_type) +
  # theme_minima
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
  theme(legend.position = c(0.837, 0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) ### PLOT SIZE - w 700 x h 500



### Splice Factor portion analysis
hepg2_k562_replicate_result <- do.call(rbind, hepg2_k562_benchmark_result_list)
k562_hepg2_replicate_result <- do.call(rbind, k562_hepg2_benchmark_result_list)

hepg2_k562_replicate_result$sig <- ifelse(hepg2_k562_replicate_result$gsea_rank <= 10, "sig", "notsig")
k562_hepg2_replicate_result$sig <- ifelse(k562_hepg2_replicate_result$gsea_rank <= 10, "sig", "notsig")

hepg2_k562_replicate_result$spl_factor <- ifelse(hepg2_k562_replicate_result$rbp_s1 %in% sf_list, "SF", "notSF")
k562_hepg2_replicate_result$spl_factor <- ifelse(k562_hepg2_replicate_result$rbp_s1 %in% sf_list, "SF", "notSF")

hepg2_k562_replicate_result$study <- "HepG2 using K562 Library"
k562_hepg2_replicate_result$study <- "K562 using HepG2 Library"

combined_replicate_result <- rbind(hepg2_k562_replicate_result, k562_hepg2_replicate_result)
combined_replicate_result$spl_type <- factor(combined_replicate_result$spl_type,
                                             levels = spl_types)

### plot for select top
count_df <- as.data.frame(table(combined_replicate_result[,c('sig', 'rbp_s1')]))
sig_rbps <- as.character(subset(count_df[count_df$sig == "sig", ], Freq > 3)$rbp_s1)
sig_rbps_sorted <- sig_rbps[order(subset(count_df[count_df$sig == "sig", ], Freq > 3)$Freq)]

heatmap_plot_df <- subset(combined_replicate_result, rbp_s1 %in% sig_rbps_sorted)
heatmap_plot_df$rbp_s1 <- factor(heatmap_plot_df$rbp_s1, levels = sig_rbps_sorted)


label_plot_df <- data.frame(rbp = sig_rbps_sorted,
                            sig = ifelse(sig_rbps_sorted %in% sf_list, "SF", "notSF"))
label_plot_df$rbp <- factor(label_plot_df$rbp, levels = sig_rbps_sorted)


ggplot(heatmap_plot_df,
       aes(x = spl_type, fill = sig, y = rbp_s1)) +
  geom_raster() +
  facet_wrap(~study) +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("sig" = "indianred2",
                               "notsig" = "grey90"),
                    labels = c("sig" = "Reproducible",
                               "notsig" = "Not reproducible")) +
  labs(fill = "Reproducible\n(Rank <= 10)",
       y = "shRNA-Perturbed RBPs",
       x = "Splicing Type") +
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
  )  +
  theme(legend.position = c(0.753, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.7), color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) ### SIZE - w 400px * h 500px




###### TEMP-DEV #####

# gather relevant data
hepg2_cdf_out <- subset(hepg2_top_rank_cdf_df, rank %in% c(1, 3, 5, 10))
k562_cdf_out <- subset(k562_top_rank_cdf_df, rank %in% c(1, 3, 5, 10))

# annotate study
hepg2_cdf_out$target <- "HepG2"
hepg2_cdf_out$library <- "K562"
k562_cdf_out$target <- "K562"
k562_cdf_out$library <- "HepG2"

out_df <- rbind(hepg2_cdf_out, k562_cdf_out)

# annotate benchmark condition
out_df$threshold <- count_threshold
out_df$count_method <- test_method


hepg2_plot_per_method[[test_method]] <<- hepg2_plot_list
k562_plot_per_method[[test_method]] <<- k562_plot_list



benchmark_df <- do.call(rbind, rank_result_df)

cowplot::plot_grid(p_hepg2 + theme(legend.position = "none"),
                   p_k562 + theme(legend.position = "bottom",
                                  legend.direction = "vertical"),
                   nrow = 1,
                   align= "hv", axis= "tb")



ggplot(subset(benchmark_df, rank == 1), aes(fill = rank, group = rank, y = count, x = target)) + geom_bar(stat = "identity", position = "dodge") + facet_grid(threshold~method)



ggplot(benchmark_df, aes(x = rank, y = count, color = method)) + geom_point() + geom_line() + facet_grid(count_method + target ~ threshold)


ggplot(benchmark_df, aes(x = factor(rank), y = count, fill = method, group = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(count_method + target ~ threshold)


# pulled PSI benchmark
ggplot(subset(benchmark_df, count_method %in% c("average", "indiv_psi") & threshold == 20), aes(x = factor(rank), y = count, fill = method, group = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(labels = c("concord_rank" = "Kendall Correlation Rank",
                                 "gsea_rank" = "GSEA Rank",
                                 "lincor_rank" = "Pearson Correlation Rank")) +
  facet_grid(target ~ count_method,
             labeller = labeller(count_method = c("average" = "Pulled PSI",
                                                  "indiv_psi" = "Mean PSI"))) +
  labs(x = "Rank of Target RBP",
       y = "Count for the Tested Cell Line",
       fill = "Correlation Methods")


# Count Benchmark
benchmark_df_count_df <- subset(benchmark_df, count_method == "average")
ggplot(benchmark_df_count_df,
       aes(x = factor(rank), y = count, fill = method, group = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(labels = c("concord_rank" = "Kendall Correlation Rank",
                                 "gsea_rank" = "GSEA Rank",
                                 "lincor_rank" = "Pearson Correlation Rank")) +
  facet_grid(target ~ threshold,
             labeller = labeller(count_method = c("average" = "Pulled PSI",
                                                  "indiv_psi" = "Mean PSI"))) +
  geom_hline(yintercept = 42, color = "indianred2", linetype = "dashed") +
  labs(x = "Rank of Target RBP",
       y = "Count for the Tested Cell Line",
       fill = "Correlation Methods")


# FDR Benchmark
benchmark_df_fdr_df <- subset(benchmark_df, count_method %in% c("average", "average.binom_fdr") & threshold == 20)
ggplot(benchmark_df_fdr_df,
       aes(x = factor(rank), y = count, fill = method, group = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(labels = c("concord_rank" = "Kendall Correlation Rank",
                                 "gsea_rank" = "GSEA Rank",
                                 "lincor_rank" = "Pearson Correlation Rank")) +
  facet_grid(target ~ count_method,
             labeller = labeller(count_method = c("average" = "Pulled PSI",
                                                  "indiv_psi" = "Mean PSI",
                                                  "average.binom_fdr" = "Pulled PSI with Binomial FDR < 0.05"))) +
  geom_hline(yintercept = 43, color = "indianred2", linetype = "dashed") +
  labs(x = "Rank of Target RBP",
       y = "Count for the Tested Cell Line",
       fill = "Correlation Methods",
       subtitle = "RC >= 20")

