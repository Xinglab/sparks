##### CODE #####
library(SPARKS)
library(ggplot2)
library(dplyr)
library(data.table)

spl_types <- c("SE", "A3SS", "A5SS")

##### ANALYSIS #####
merged_kd_library <- readRDS('/Users/harryyang/research/Xing_Lab/sparks/data/ENCODE_library.KD_and_KO.count_20.rds')

# determine replicate RBP experiments
shrna_study_list <- names(merged_kd_library$SE)[grep("shRNA", names(merged_kd_library$SE))]

# generate study -> RBP table
shrna_rbp_info_df <- data.frame(study = shrna_study_list,
                                rbp = unlist(lapply(shrna_study_list, function(x) strsplit(x, "_")[[1]][2])))
shrna_rbp_num_exp <- table(shrna_rbp_info_df$rbp)
shrna_rbp_kd2_list <- names(shrna_rbp_num_exp)[shrna_rbp_num_exp == 2]

# determine replicate RBP experiments
crispr_study_list <- names(merged_kd_library$SE)[grep("CRISPR", names(merged_kd_library$SE))]

# generate study -> RBP table
crispr_rbp_info_df <- data.frame(study = crispr_study_list,
                                 rbp = unlist(lapply(crispr_study_list, function(x) strsplit(x, "_")[[1]][2])))
crispr_rbp_num_exp <- table(crispr_rbp_info_df$rbp)
crispr_rbp_kd2_list <- names(crispr_rbp_num_exp)[crispr_rbp_num_exp == 2]

##### SUMMARY TABLE #####
num_shrna_k562 <- length(grep(shrna_study_list, pattern = "K562"))
num_shrna_hepg2 <- length(grep(shrna_study_list, pattern = "HepG2"))
num_crispr_k562 <- length(grep(crispr_study_list, pattern = "K562"))
num_crispr_hepg2 <- length(grep(crispr_study_list, pattern = "HepG2"))

num_shrna_common <- length(shrna_rbp_kd2_list)
num_crispr_common <- length(crispr_rbp_kd2_list)

rbps_in_shrna_crispr <- intersect(unique(shrna_rbp_info_df$rbp), unique(crispr_rbp_info_df$rbp))
num_rbp_in_shrna_crispr <- length(rbps_in_shrna_crispr)

##### DEV ####
benchmark_result_df <- data.table::fread("/Users/harryyang/Documents/Research/Xing_Lab/sparks/data/ENCODE.pairwise_all.SE.df.txt")

# annotate cell lin and RBP for easier processing
benchmark_result_df$study_s1 <- unlist(lapply(benchmark_result_df$S1, function(x) strsplit(x, "_")[[1]][1]))
benchmark_result_df$study_s2 <- unlist(lapply(benchmark_result_df$S2, function(x) strsplit(x, "_")[[1]][1]))
benchmark_result_df$rbp_s1 <- unlist(lapply(benchmark_result_df$S1, function(x) strsplit(x, "_")[[1]][2]))
benchmark_result_df$rbp_s2 <- unlist(lapply(benchmark_result_df$S2, function(x) strsplit(x, "_")[[1]][2]))
benchmark_result_df$method_s1 <- unlist(lapply(benchmark_result_df$S1, function(x) strsplit(x, "_")[[1]][3]))
benchmark_result_df$method_s2 <- unlist(lapply(benchmark_result_df$S2, function(x) strsplit(x, "_")[[1]][3]))



replicate_subset <- subset(benchmark_result_df, rbp_s1 %in% rbps_in_shrna_crispr & rbp_s2 %in% rbps_in_shrna_crispr)

## generate summary for number of category items
# category 1 - same RBP perturbation in same cell line but different method
# category 2 - same RBP perturbation in different cell lien and different method

View(table(unique(replicate_subset[, c('study_s1','rbp_s1','method_s1')])))

experiment_table <- as.data.frame(table(unique(replicate_subset[, c('study_s1','rbp_s1','method_s1')])))

combo_guide_list <- list()

experiment_table %>%
  filter(Freq == 1) %>%  # keep the experiments with any values
  group_by(rbp_s1,
           .drop = F) %>%
  group_map( ~ {  # apply the function to generate combinations
    shrna <- subset(.x, method_s1 == "shRNA") %>% arrange(study_s1)
    crispr <- subset(.x, method_s1 == "CRISPR") %>% arrange(study_s1)
    # print(shrna)
    # print(crispr)

    # generate combinations
    combo_table <- data.table::CJ(seq(dim(shrna)[1]), seq(dim(crispr)[1]))
    dummy <- apply(combo_table, 1, function(a){
      combo_guide_entry <- t(cbind(c(shrna[a[1], c('rbp_s1', 'study_s1', 'method_s1'), drop = F], # since the rbp is the same, we are extracting them here
                                                   crispr[a[2], c('study_s1', 'method_s1'), drop = F])))
      # rename the columns
      # print(combo_guide_entry)

      combo_guide_entry <- lapply(combo_guide_entry, as.character)
      names(combo_guide_entry) <- c("rbp", "study_s1", "method_s1", "study_s2", "method_s2")

      combo_guide_list[[length(combo_guide_list) + 1]] <<- combo_guide_entry
    })
    return()
  },
  .keep = T)

combo_guide_table_temp <- data.frame(do.call(rbind, combo_guide_list))  # TODO - this returns each column as list - not sure why
combo_guide_table <- do.call(rbind, apply(combo_guide_table_temp, 1, function(x) data.frame(x)))  # TODO - refactor the above code

# subset the same cell line vs. other cell line
combo_same_cell <- subset(combo_guide_table, study_s1 == study_s2)
combo_diff_cell <- subset(combo_guide_table, study_s1 != study_s2)


# calculate rank
perform_cross_method_benchmark_calculation <- function(entry){
  # test one direction
  test_df <- subset(replicate_subset,
                    rbp_s1 == entry['rbp'] & study_s1 == entry['study_s1'] & method_s1 == entry['method_s1'] &  # select the entries for the study
                      study_s2 == entry['study_s2'] & method_s2 == entry['method_s2'])
  result_df <- calculate_benchmark_rank(test_df)
  benchmark_result_df <- extract_benchmark_replicate_result(result_df)
  combined_result_df_one <- cbind(benchmark_result_df, t(entry))

  # test other direction
  # flip the entry
  entry_rev <- entry
  entry_rev['study_s1'] <- entry['study_s2']
  entry_rev['study_s2'] <- entry['study_s1']
  entry_rev['method_s1'] <- entry['method_s2']
  entry_rev['method_s2'] <- entry['method_s1']

  test_df_rev <- subset(replicate_subset,
                    rbp_s1 == entry_rev['rbp'] & study_s1 == entry_rev['study_s1'] & method_s1 == entry_rev['method_s1'] &  # select the entries for the study
                      study_s2 == entry_rev['study_s2'] & method_s2 == entry_rev['method_s2'])
  result_df_rev <- calculate_benchmark_rank(test_df_rev)
  benchmark_result_df_rev <- extract_benchmark_replicate_result(result_df_rev)

  combined_result_df_rev <- cbind(benchmark_result_df_rev, t(entry_rev))

  # combine the result from both test
  combined_result_df <- rbind(combined_result_df_one, combined_result_df_rev)
  return(combined_result_df)
}


same_cell_result <- do.call(rbind, apply(combo_same_cell, 1, function(entry) perform_cross_method_benchmark_calculation(entry)))
diff_cell_result <- do.call(rbind, apply(combo_diff_cell, 1, function(entry) perform_cross_method_benchmark_calculation(entry)))

same_cell_cdf_crispr_to_shrna <- perform_cross_method_cdf_calculation(subset(same_cell_result, method_s1 == "CRISPR"))
same_cell_cdf_shrna_to_crispr <- perform_cross_method_cdf_calculation(subset(same_cell_result, method_s1 == "shRNA"))
diff_cell_cdf_crispr_to_shrna <- perform_cross_method_cdf_calculation(subset(diff_cell_result, method_s1 == "CRISPR"))
diff_cell_cdf_shrna_to_crispr <- perform_cross_method_cdf_calculation(subset(diff_cell_result, method_s1 == "shRNA"))

# add labels
same_cell_cdf_crispr_to_shrna$anno <- "crispr"
same_cell_cdf_shrna_to_crispr$anno <- "shrna"
diff_cell_cdf_crispr_to_shrna$anno <- "crispr"
diff_cell_cdf_shrna_to_crispr$anno <- "shrna"

same_cell_cdf_combined <- rbind(same_cell_cdf_crispr_to_shrna,
                                same_cell_cdf_shrna_to_crispr)
diff_cell_cdf_combined <- rbind(diff_cell_cdf_crispr_to_shrna,
                                diff_cell_cdf_shrna_to_crispr)

same_cell_cdf_combined$diff <- 'same'
diff_cell_cdf_combined$diff <- 'diff'

cross_cdf_combined <- rbind(same_cell_cdf_combined,
                            diff_cell_cdf_combined)
# generate plot
p_benchmark <- ggplot(cross_cdf_combined,
       aes(x = rank,
           y = count,
           color = method)) +
  geom_point() +
  geom_line(aes(linetype = anno)) +
  scale_color_discrete(labels = new_test_names) +
  scale_x_continuous(breaks = seq(-1, 10),
                     limits = c(1, 10)) +
  scale_y_continuous(limits = c(0, 40), breaks = c(seq(0, 5) * 8)) +  # scale for shRNA KD
  # scale_y_continuous(limits = c(0, 20), breaks = c(seq(0, 5) * 4)) +  # scale for same cell - max = 24
  scale_linetype_discrete(labels = c("shrna" = 'shRNA using CRISPR library',
                                     'crispr' = "CRISPR using shRNA library")) +
  labs(x = "Rank of the target RBP from each method",
       y = "Cumulative Sum of RBP KD Experiments",
       color = "Method",
       linetype = "Benchmark Strategy") +
  facet_grid(~diff,
             labeller = labeller(diff = c("diff" = "Different cell line (n = 40)",
                                          "same" = "Same cell line (n = 24)"))) +
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
  )

cowplot::plot_grid(p_benchmark +
                     theme(legend.position = "none"),
                   cowplot::get_legend(p_benchmark +
                     theme(legend.background = element_rect(fill = "white", color = "black"),
                           legend.text = element_text(size = 8),
                           legend.title = element_text(size = 8),
                           legend.position = "right")),
                   rel_widths = c(3, 1))



perform_cross_method_cdf_calculation <- function(result){
  aa <- as.data.frame(result)
  # print(aa)
  cdf_df <- do.call(rbind, lapply(c("lincor_rank", "concord_rank", "gsea_rank"), function(x) {

    cumsum_df <- as.data.frame(unlist(lapply(seq(unique(aa$rbp)), function(y)
      return(sum(aa[, x] <= y ))
    )))

    colnames(cumsum_df) <- c("cum_count")
    res_df <- data.frame(rank = as.numeric(rownames(cumsum_df)),
                         count = cumsum_df$cum_count,
                         method = x)

    return(res_df)
  }))
  return(cdf_df)
}



test_df <- subset(replicate_subset, S1 == "HepG2_NFX1_CRISPR" & method_s2 == "shRNA" & study_s2 == "K562")
View(calculate_benchmark_rank(test_df))



