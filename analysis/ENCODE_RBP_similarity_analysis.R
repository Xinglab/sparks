##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(pbmcapply)
library(igraph)
library(ggVennDiagram)
library(dplyr)

##### FUNCTIONS #####
# calculate rank
calculate_benchmark_rank <- function(result_df){

  # process each target study result separately
  target_study_list <- unique(result_df$S1)

  rank_anno_result_df <- do.call(rbind, lapply(target_study_list, function(target_study){
    target_result_df <- subset(result_df, S1 == target_study)


    # calculate rank
    target_result_df$lincor_rank <- rank(-target_result_df$score_abs)
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
  rep_plot_df <- rep_result[, c("lincor_rank", "concord_rank", "gsea_rank", "gsea_score", "rbp_s1")]
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


# extract merged peaks
extract_merged_peaks <- function(test_clip_exp){
  # TODO - fix this hard code below - may not work for A3SS/A5SS
  clip_region_list <- names(clip_library)
  clip_merged_df <- do.call(cbind, lapply(clip_region_list, function(clip_region){
    clip_region_df <- clip_library[[clip_region]][, test_clip_exp, drop = F]
    return(clip_region_df)
  }))
  clip_merged <- ifelse(rowSums(clip_merged_df == "peak") > 0, "peak", "no_peak")
  # clip_merged <- factor(clip_merged, levels = c("peak", "no_peak"))
  names(clip_merged) <- rownames(clip_merged_df)
  return(clip_merged)
}


# calculate CLIP jaccard index
calculate_jaccard_index_for_clip <- function(clip_list, target_cell){
  clip_jaccard_index <- do.call(rbind, apply(t(combn(clip_list, 2)), 1, function(x){
    print(x)

    exp_1 <- x[1]
    exp_2 <- x[2]
    spl_exp_1 <- sprintf("%s_%s_shRNA", target_cell, exp_1)
    spl_exp_2 <- sprintf("%s_%s_shRNA", target_cell, exp_2)
    # extract events
    exp_1_events <- extract_GSEA_significant_events(kd_library[[spl_exp_1]])
    exp_2_events <- extract_GSEA_significant_events(kd_library[[spl_exp_2]])

    exp_1_events_pos <- exp_1_events$positive
    exp_1_events_neg <- exp_1_events$negative
    exp_2_events_pos <- exp_2_events$positive
    exp_2_events_neg <- exp_2_events$negative

    # extract concordant events
    concordant_events <- union(intersect(exp_1_events_pos, exp_2_events_pos),
                               intersect(exp_1_events_neg, exp_2_events_neg))

    ## extract clip data
    # merge the clip peaks for all regions
    test_clip_1 <- sprintf("%s_%s", target_cell, exp_1)
    test_clip_2 <- sprintf("%s_%s", target_cell, exp_2)

    clip_1_merged <- extract_merged_peaks(test_clip_1)
    clip_2_merged <- extract_merged_peaks(test_clip_2)

    # combine the clip peak results from the two results
    merged_clip_df <- as.data.frame(cbind(clip_1_merged, clip_2_merged))

    merged_clip_df_select <- merged_clip_df[concordant_events, ]
    # merged_clip_df_select <- merged_clip_df


    merged_clip_df_clean <- merged_clip_df_select[rowSums(is.na(merged_clip_df_select)) == 0, ]

    # calculate the counts
    # table(merged_clip_df_clean) - this approach does not work if one CLIP exp is totally empty
    count_no_no <- dim(subset(merged_clip_df_clean, clip_1_merged == "no_peak" & clip_2_merged == "no_peak"))[1]
    count_no_yes <- dim(subset(merged_clip_df_clean, clip_1_merged == "no_peak" & clip_2_merged == "peak"))[1]
    count_yes_no <- dim(subset(merged_clip_df_clean, clip_1_merged == "peak" & clip_2_merged == "no_peak"))[1]
    count_yes_yes <- dim(subset(merged_clip_df_clean, clip_1_merged == "peak" & clip_2_merged == "peak"))[1]

    # calculate jaccard index
    jaccard_index <- count_yes_yes / (count_yes_yes + count_yes_no + count_no_yes)

    # make it into a table
    result_df <- data.frame(rbp_1 = exp_1,
                            rbp_2 = exp_2,
                            jaccard = jaccard_index,
                            no_no = count_no_no,
                            no_yes = count_no_yes,
                            yes_no = count_yes_no,
                            yes_yes = count_yes_yes,
                            distance = distances(combined_graph, exp_1, exp_2)[1],
                            cell_line = target_cell)
    return(result_df)
  }))
  return(clip_jaccard_index)
}


# extract adjacency matrix
generate_adjacency_matrix_from_benchmark <- function(result_df_anno, score_threshold = 0.6){
  # extract top results
  top_results <- subset(result_df_anno, abs(gsea_score) >= score_threshold)

  # extract top experiments for matrix construction
  interest_exp_list <- unique(c(top_results$S1, top_results$S2))

  # extract corresponding rbp list
  rbp_list <- sort(unique(c(top_results$rbp_s1, top_results$rbp_s2)))

  ### construct adjacency matrix
  # make empty adjacency matrix
  adj_matrix <- matrix(0, nrow = length(rbp_list), ncol = length(rbp_list))
  colnames(adj_matrix) <- rbp_list
  rownames(adj_matrix) <- rbp_list

  # fill in the values

  dummy <- apply(top_results, 1, function(entry){
    adj_matrix[entry["rbp_s1"], entry["rbp_s2"]] <<- adj_matrix[entry["rbp_s1"], entry["rbp_s2"]] + 1
    return()
  })

  return(adj_matrix)
}


generate_score_matrix_from_benchmark <- function(result_df_anno, rbp_list, score_threshold = 0.6){
  # extract top results
  top_results <- subset(result_df_anno, gsea_score >= score_threshold)

  # extract top experiments for matrix construction
  interest_exp_list <- unique(c(top_results$S1, top_results$S2))

  # extract corresponding rbp list
  # rbp_list <- sort(unique(c(top_results$rbp_s1, top_results$rbp_s2)))

  ### construct adjacency matrix
  # make empty adjacency matrix
  score_matrix <- matrix(0, nrow = length(rbp_list), ncol = length(rbp_list))
  colnames(score_matrix) <- rbp_list
  rownames(score_matrix) <- rbp_list

  # fill in the values
  all_results <- subset(result_df_anno, S1 %in% interest_exp_list & S2 %in% interest_exp_list)
  dummy <- apply(all_results, 1, function(entry){
    score_matrix[entry["rbp_s1"], entry["rbp_s2"]] <<- as.numeric(entry['gsea_score'])
    return()
  })

  return(score_matrix)
}


calculate_network_binary_df <- function(kd_library, test_cell, test_rbps, test_method = "shRNA"){
  test_list <- unlist(lapply(test_rbps, function(x) sprintf("%s_%s_%s", test_cell, x, test_method)))

  test_total_events <- unique(unlist(lapply(test_list, function(x) unlist(extract_GSEA_significant_events(kd_library[[x]])))))

  # gather directional data
  test_network_binary_df <- as.data.frame(do.call(cbind, lapply(test_list, function(x){
    sig_events <- extract_GSEA_significant_events(kd_library[[x]])
    sig_binary_list <- ifelse(test_total_events %in% sig_events$positive, 1,
                              ifelse(test_total_events %in% sig_events$negative, -1, 0))

  })))

  # polish up the binary df
  rownames(test_network_binary_df) <- test_total_events
  colnames(test_network_binary_df) <- test_rbps

  return(test_network_binary_df)
}


## generate count bar plot
generate_event_count_plot <- function(input_network_df, cell_line){
  event_count_df <- data.frame(table(reshape2::melt(input_network_df)))

  p <- ggplot(subset(event_count_df, value != 0), aes(x = variable, fill = value, y = Freq)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c('-1' = 'lightcoral',
                                 '1' = 'skyblue1'),
                      labels = c('-1' = "More Included in KD",
                                 '1' = "More Excluded in KD")) +
    theme_bw() +
    scale_y_continuous(limits = c(0, ceiling(max(colSums(abs(input_network_df)))/250) * 250),
                       expand = expansion(0, 0)) +
    labs(fill = "Splicing Direction",
         y = "Number of SE events",
         x = sprintf("Perturbed RBP in %s", cell_line)) +
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    # theme_minima
    theme(axis.text.y   = element_text(size=10),
          axis.text.x   = element_text(size=10),
          axis.title.y  = element_text(size=12),
          axis.title.x  = element_text(size=12),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.line = element_line(colour = "black"),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1)
    )
  return(p)
}


import_benchmark_df_list <- function(spl_type){
  # read in the benchmark result file
  input_file <- sprintf('/Users/harryyang/Documents/Research/Xing_Lab/sparks/data/benchmark_result.shRNA.%s.df.txt', spl_type)
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

  # gather results for export
  result_df_list <- list()
  result_df_list[['HepG2_K562']] <- hepg2_k562_result_df_anno
  result_df_list[['K562_HepG2']] <- k562_hepg2_result_df_anno
  result_df_list[['HepG2_HepG2']] <- hepg2_hepg2_result_df
  result_df_list[['K562_K562']] <- k562_k562_result_df
  return(result_df_list)
}


compute_rank_df_list <- function(benchmark_df_list){
  # load them
  hepg2_k562_result_df_anno <- benchmark_df_list$HepG2_K562
  k562_hepg2_result_df_anno <- benchmark_df_list$K562_HepG2

  ## extract the rank for replicate for downstream analysis
  hepg2_k562_rank_df <- extract_benchmark_replicate_result(hepg2_k562_result_df_anno)
  k562_hepg2_rank_df <- extract_benchmark_replicate_result(k562_hepg2_result_df_anno)

  # annotate the study
  hepg2_k562_rank_df$study <- "HepG2 using K562 Library"
  k562_hepg2_rank_df$study <- "K562 using HepG2 Library"

  # gather results for export
  result_df_list <- list()
  result_df_list[['HepG2_K562']] <- hepg2_k562_rank_df
  result_df_list[['K562_HepG2']] <- k562_hepg2_rank_df

  return(result_df_list)
}


extract_reproducible_rbps <- function(rank_df_list){
  # load the elements
  hepg2_k562_rank_df <- rank_df_list$HepG2_K562
  k562_hepg2_rank_df <- rank_df_list$K562_HepG2

  reproducible_rbps <-intersect(subset(k562_hepg2_rank_df, gsea_rank <= 10)$rbp,
                                 subset(hepg2_k562_rank_df, gsea_rank <= 10)$rbp)
  return(reproducible_rbps)
}


## Generate enrichment score distribution plot
generate_enrichment_score_distribution_plot <- function(rank_df_list,
                                                        reproducible_rbps,
                                                        score_threshold = 0.6){

  # merge the score df
  cross_cell_rank_df <- do.call(rbind, rank_df_list)
  cross_cell_rank_df$reprod <- factor(ifelse(cross_cell_rank_df$rbp_s1 %in% reproducible_rbps,
                                             "Reproducible",
                                             "Not Reproducible"),
                                      levels = c("Reproducible",  # define column orders
                                                 "Not Reproducible"))

  p <- ggplot(cross_cell_rank_df, aes(x = reprod,
                                      y = gsea_score,
                                      fill = reprod)) +
    geom_hline(yintercept = score_threshold,  # add threshold line
               color = "indianred3",
               linetype = "dashed") +
    geom_boxplot(width = 0.5,
                 alpha = 0.8,
                 outlier.shape = NA) +
    geom_jitter(width = 0.25,
                alpha = 0.5) +
    scale_fill_manual(values = c("Reproducible" = "dodgerblue",
                                 "Not Reproducible" = "paleturquoise")) +
    scale_x_discrete(labels = c("Reproducible" = sprintf("Reproducible\n(n = %s)",
                                                         length(reproducible_rbps)),
                                "Not Reproducible" = sprintf("Not Reproducible\n(n = %s)",
                                                             sum(!unique(cross_cell_rank_df$rbp_s1) %in% reproducible_rbps)))) +
    scale_y_continuous(limits = c(floor(min(cross_cell_rank_df$gsea_score)/0.3),  # increment of 0.3 in PSI
                                  ceiling(max(cross_cell_rank_df$gsea_score)/0.3)) * 0.3,
                       breaks = c(seq(floor(min(cross_cell_rank_df$gsea_score)/0.3),
                                      ceiling(max(cross_cell_rank_df$gsea_score)/0.3))) * 0.3) +
    facet_wrap(~ study) +
    theme_bw() +
    labs(x = "Reproducible between two cell lines (rank <= 10)",
         y = "Enrichment Score") +
    theme(legend.position = "none") +
    theme(axis.text.y   = element_text(size = 8),
          axis.text.x   = element_text(size = 8),
          axis.title.y  = element_text(size = 10),
          axis.title.x  = element_text(size = 10),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.line = element_line(colour = "black"),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1)
    )
  return(p)
}



# calculate benchmark rank for cross cell line, while removing duplicate runs - this will have highest rank
extract_within_cell_line_rank_df <- function(benchmark_df_list, reproducible_rbps){
  # load them
  hepg2_hepg2_result_df <- benchmark_df_list$HepG2_HepG2
  k562_k562_result_df <- benchmark_df_list$K562_K562

  # clean up the data
  hepg2_result_df_anno <- calculate_benchmark_rank(subset(hepg2_hepg2_result_df, rbp_s1 != rbp_s2 & rbp_s1 %in% reproducible_rbps & rbp_s2 %in% reproducible_rbps))
  k562_result_df_anno <- calculate_benchmark_rank(subset(k562_k562_result_df, rbp_s1 != rbp_s2 & rbp_s1 %in% reproducible_rbps & rbp_s2 %in% reproducible_rbps))
  # DEV CODE - apply to all
  # hepg2_result_df_anno <- calculate_benchmark_rank(subset(hepg2_hepg2_result_df, rbp_s1 != rbp_s2))
  # k562_result_df_anno <- calculate_benchmark_rank(subset(k562_k562_result_df, rbp_s1 != rbp_s2))

  # gather results for export
  result_df_list <- list()
  result_df_list[['HepG2_HepG2']] <- hepg2_result_df_anno
  result_df_list[['K562_K562']] <- k562_result_df_anno

  return(result_df_list)
}


calculate_combined_avg_enrichment_score <- function(win_cl_rank_df, score_threshold = 0.6){
  # import data
  hepg2_result_df_anno <- win_cl_rank_df$HepG2_HepG2
  k562_result_df_anno <- win_cl_rank_df$K562_K562

  # extract RBP list for fitted list
  rbp_list <- sort(unique(c(hepg2_result_df_anno$rbp_s1,
                       hepg2_result_df_anno$rbp_s2,
                       k562_result_df_anno$rbp_s1,
                       k562_result_df_anno$rbp_s2)))

  # calculate score
  k562_score_matrix <- generate_score_matrix_from_benchmark(k562_result_df_anno,
                                                            rbp_list,
                                                            score_threshold = score_threshold)

  hepg2_score_matrix <- generate_score_matrix_from_benchmark(hepg2_result_df_anno,
                                                             rbp_list,
                                                             score_threshold = score_threshold)

  # take average
  hepg2_score_avg <- (hepg2_score_matrix + t(hepg2_score_matrix)) / 2
  hepg2_score_avg[hepg2_score_avg < score_threshold] <- 0

  k562_score_avg <- (k562_score_matrix + t(k562_score_matrix)) / 2
  k562_score_avg[k562_score_avg < score_threshold] <- 0

  common_rbps <- sort(intersect(rownames(k562_score_avg), rownames(hepg2_score_avg)))

  combined_avg <- (k562_score_avg[common_rbps, common_rbps] + hepg2_score_avg[common_rbps, common_rbps]) / 2
  combined_avg[combined_avg < score_threshold] <- 0
  return(combined_avg)
}


##### DATA PREPARATION #####
count_threshold <- 20
## process cluster-run benchmark data
# read in benchmark result data

# import benchmark results
benchmark_df_list_se <- import_benchmark_df_list("SE")
benchmark_df_list_a5ss <- import_benchmark_df_list("A5SS")
benchmark_df_list_a3ss <- import_benchmark_df_list("A3SS")

# calcualte rank
rank_df_list_se <- compute_rank_df_list(benchmark_df_list_se)
rank_df_list_a5ss <- compute_rank_df_list(benchmark_df_list_a5ss)
rank_df_list_a3ss <- compute_rank_df_list(benchmark_df_list_a3ss)

# find reproducible RBPs in each
reproducible_rbps_se <- extract_reproducible_rbps(rank_df_list_se)
reproducible_rbps_a3ss <- extract_reproducible_rbps(rank_df_list_a3ss)
reproducible_rbps_a5ss <- extract_reproducible_rbps(rank_df_list_a5ss)

# generate enrichment score distribution plot
score_threshold = 0.6
generate_enrichment_score_distribution_plot(rank_df_list_se, reproducible_rbps_se, score_threshold)
generate_enrichment_score_distribution_plot(rank_df_list_a5ss, reproducible_rbps_a5ss, score_threshold)
generate_enrichment_score_distribution_plot(rank_df_list_a3ss, reproducible_rbps_a3ss, score_threshold)

# union of significant RBPs
reproducible_rbps_combined <- unique(c(reproducible_rbps_se,
                                       reproducible_rbps_a3ss,
                                       reproducible_rbps_a5ss))

## merge splicing type results
# calculate within cell line rank df
win_cl_rank_df_se <- extract_within_cell_line_rank_df(benchmark_df_list_se, reproducible_rbps_combined)
win_cl_rank_df_a5ss <- extract_within_cell_line_rank_df(benchmark_df_list_a5ss, reproducible_rbps_combined)
win_cl_rank_df_a3ss <- extract_within_cell_line_rank_df(benchmark_df_list_a3ss, reproducible_rbps_combined)

# calculate average scores between pairwise relationship for HepG2 and K562 within cell line results
combined_avg_se <- calculate_combined_avg_enrichment_score(win_cl_rank_df_se, score_threshold)
combined_avg_a5ss <- calculate_combined_avg_enrichment_score(win_cl_rank_df_a5ss, score_threshold)
combined_avg_a3ss <- calculate_combined_avg_enrichment_score(win_cl_rank_df_a3ss, score_threshold)

# combine and average the pairwise scores in all splicing types
combined_all <- (combined_avg_se + combined_avg_a5ss + combined_avg_a3ss) / 3

# flatten the scores that are below the threshold
combined_all[combined_all < score_threshold] <- 0

# remove empty rows for better visualization
empty_rbps <- names(rowSums(combined_all))[rowSums(combined_all) == 0]
combined_filtered <- combined_all[!(rownames(combined_all) %in% empty_rbps),
                                  !(colnames(combined_all) %in% empty_rbps)]

# generate heatmap
pheatmap::pheatmap(combined_filtered,
                   color = colorRampPalette(c("white", "gold", "red"))(100),
                   clustering_distance_cols = "minkowski",
                   clustering_distance_rows = "minkowski",
                   # clustering_distance_cols = "manhattan",
                   # clustering_distance_rows = "manhattan",
                   # clustering_distance_cols = "binary",
                   # clustering_distance_rows = "binary",
                   clustering_method = "ward.D2",
                   # clustering_method = "ward.D",
                   # clustering_method = "mcquitty",

                   cutree_rows = 8,
                   cutree_cols = 8,
                   # fontsize_row = 6,
                   # fontsize_col = 6,
                   # fontsize_number = 6
                   fontsize = 8
)


##### SF3 Network Analysis
combined_adj_matrix <- combined_all
combined_adj_matrix[is.na(combined_adj_matrix)] <- 0
combined_adj_matrix[combined_adj_matrix > 0] <- 1

# select_rbps <- colnames(combined_adj_matrix)[colSums(combined_adj_matrix[c("SF3A3", "U2AF1"), ]) > 0]  # this is okay as the adj matrix is symmetrical
select_rbps <- colnames(combined_adj_matrix)[colSums(combined_adj_matrix[c("SF3A3", "SF3B1", "SF3B4", "U2AF1", "U2AF2"), ]) > 0]  # this is okay as the adj matrix is symmetrical

select_adj_matrix <- combined_adj_matrix[select_rbps, select_rbps]
select_graph <- graph_from_adjacency_matrix(select_adj_matrix, mode = "min", weighted = NULL)

# manually set color
vertex_colors <- as.vector(rep("black", length(rownames(select_adj_matrix))))
names(vertex_colors) <- rownames(select_adj_matrix)

vertex_colors[grep("SF3", names(vertex_colors))] <- "green4"
vertex_colors[grep("U2", names(vertex_colors))] <- "goldenrod3"
plot.igraph(select_graph,
            # layout = layout.fruchterman.reingold(select_graph),
            # layout = layout_nicely(select_graph),
            layout = layout_with_kk(select_graph),
            #
            # vertex.label.cex = 0.75,
            # vertex.size = 30,
            # vertex.color = 'dodgerblue1',
            vertex.label.color = vertex_colors,
            # vertex.label.family="Helvetica",
            vertex.label.family = "Arial",
            edge.arrow.size = 0.5,
            vertex.label.cex = 0.75,
            vertex.label.font = 2,
            vertex.shape = "circle",
            vertex.size = 0.1,
            vertex.label.color = "black",
            edge.width = 1)


# select_rbps <- colnames(combined_adj_matrix)[colSums(combined_adj_matrix[c("SF3A3", "U2AF1"), ]) > 0]  # this is okay as the adj matrix is symmetrical
select_rbps <- colnames(combined_adj_matrix)[colSums(combined_adj_matrix[c("MAGOH", "KHSRP", "MATR3", "NELFE", "PTBP1", "PCBP1", "SRSF1"), ]) > 0]  # this is okay as the adj matrix is symmetrical

# select_rbps <- c("MAGOH", "MATR3", "NELFE", "PTBP1", "SRSF7", "PCBP1", "SRSF1")
select_adj_matrix <- combined_adj_matrix[select_rbps, select_rbps]
select_graph <- graph_from_adjacency_matrix(select_adj_matrix, mode = "min", weighted = NULL)

# manually set color
vertex_colors <- as.vector(rep("black", length(rownames(select_adj_matrix))))
names(vertex_colors) <- rownames(select_adj_matrix)

vertex_colors[c("MAGOH", "MATR3", "NELFE", "PTBP1", "KHSRP", "PCBP1", "SRSF1")] <- "dodgerblue1"
# vertex_colors[c("DDX21", "HNRNPK", "KHSRP", "FXR1", "EIF4A3")] <- "black"

plot.igraph(select_graph,
            # layout = layout.fruchterman.reingold(select_graph),
            # layout = layout_nicely(select_graph),
            # layout = layout_with_kk(select_graph),
            layout = layout_in_circle(select_graph),

            #
            # vertex.label.cex = 0.75,
            # vertex.size = 30,
            # vertex.color = 'dodgerblue1',
            vertex.label.color = vertex_colors,
            # vertex.label.family="Helvetica",
            vertex.label.family = "Arial",
            edge.arrow.size = 0.5,
            vertex.label.cex = 0.75,
            vertex.label.font = 2,
            vertex.shape = "circle",
            vertex.size = 0.1,
            vertex.label.color = "black",
            edge.width = 1)
