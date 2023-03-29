##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(dplyr)



kd_library <- readRDS("/Users/harryyang/transfer/Clean_minfiltered.SE.library.rds")

##### ANALYSIS #####
# read in K562 S34F data
k562_sparks_file <- '/Users/harryyang/transfer/SF_Mutation/K562_U2AF1_S34F.30.SPARKS.rds'
k562_sparks <- readRDS(k562_sparks_file)
k562_mats <- import_SPARKS_MATS_for_analysis(k562_sparks, "SE")

k562_result <- k562_sparks@SPARKS_analysis_result$SE

add_plot_title(generate_enrichment_barplot(k562_result,
                                           bar_color = "violetred1",
                                           num_plot = 10,
                                           select_genes = c("U2AF1", "PTBP1")),
               "K562 U2AF1 S34F",
               title_color = "violetred1") # 350 x 400

# read in HEK293 PTBP1 KD  data
hek_kd_sparks_file <- '/Users/harryyang/transfer/HEK/HEK293_PTBP1_KD.30.SPARKS.rds'
hek_kd_sparks <- readRDS(hek_kd_sparks_file)
hek_kd_mats <- import_SPARKS_MATS_for_analysis(hek_kd_sparks, "SE")

hek_kd_result <- hek_kd_sparks@SPARKS_analysis_result$SE

add_plot_title(generate_enrichment_barplot(hek_kd_result,
                                           bar_color = "darksalmon",
                                           num_plot = 10,
                                           select_genes = c("U2AF1", "PTBP1")),
               "HEK293 PTBP1 siRNA",
               title_color = "darksalmon") # 350 x 400

# read in HEK E9KO data
hek_e9ko_sparks_file <- '/Users/harryyang/transfer/HEK/HEK293_PTBP1_delE9rescue.30.SPARKS.rds'
hek_e9ko_sparks <- readRDS(hek_e9ko_sparks_file)
hek_e9ko_mats <- import_SPARKS_MATS_for_analysis(hek_e9ko_sparks, "SE")

hek_e9ko_result <- hek_e9ko_sparks@SPARKS_analysis_result$SE

add_plot_title(generate_enrichment_barplot(hek_e9ko_result,
                                           bar_color = "darksalmon",
                                           num_plot = 10,
                                           select_genes = c("U2AF1", "PTBP1")),
               "K562 U2AF1 S34F",
               title_color = "darksalmon") # 350 x 400

# read in K562 S34F data
hek_e9ki_sparks_file <- '/Users/harryyang/transfer/HEK/HEK293_PTBP1_FLrescue.30.SPARKS.rds'
hek_e9ki_sparks <- readRDS(hek_e9ki_sparks_file)
hek_e9ki_mats <- import_SPARKS_MATS_for_analysis(hek_e9ki_sparks, "SE")

hek_e9ki_result <- hek_e9ki_sparks@SPARKS_analysis_result$SE

add_plot_title(generate_enrichment_barplot(hek_e9ki_result,
                                           bar_color = "darksalmon",
                                           num_plot = 10,
                                           select_genes = c("U2AF1", "PTBP1")),
               "K562 U2AF1 S34F",
               title_color = "darksalmon") # 350 x 400


# read in K562 S34F data
hek_e9koki_sparks_file <- '/Users/harryyang/transfer/HEK/HEK293_PTBP1_delE9.30.SPARKS.rds'
hek_e9koki_sparks <- readRDS(hek_e9koki_sparks_file)
hek_e9koki_mats <- import_SPARKS_MATS_for_analysis(hek_e9koki_sparks, "SE")

hek_e9koki_result <- hek_e9koki_sparks@SPARKS_analysis_result$SE

add_plot_title(generate_enrichment_barplot(hek_e9koki_result,
                                           bar_color = "darksalmon",
                                           num_plot = 10,
                                           select_genes = c("U2AF1", "PTBP1")),
               "HEK293 E9KO vs. E9KI",
               title_color = "darksalmon") # 350 x 400

### HEK ANALYSIS ###
generate_RBP_KD_correlation_scatter_plot(hek_kd_mats,
                                         hek_e9ki_mats,
                                         "HEK293 PTBP1 KD vs. Control", "HEK293 PTBP1 KD vs. E9 KI Rescue",
                                         slope = T)
generate_RBP_KD_correlation_scatter_plot(hek_kd_mats,
                                         hek_e9ko_mats,
                                         "HEK293 PTBP1 KD vs. Control", "HEK293 PTBP1 KD vs. E9 KO Rescue",
                                         slope = T)
generate_RBP_KD_correlation_scatter_plot(hek_e9ki_mats,
                                         hek_e9ko_mats,
                                         "HEK293 PTBP1 KD vs. E9 KI Rescue", "HEK293 PTBP1 KD vs. E9 KO Rescue",
                                         slope = T)
generate_RBP_KD_correlation_scatter_plot(hek_kd_mats,
                                         hek_e9koki_mats,
                                         "HEK293 PTBP1 KD vs. Control", "HEK293 PTBP1 E9 KI vs. E9 KO Rescue",
                                         slope = T)

generate_RBP_KD_correlation_scatter_plot(kd_library$HepG2_PTBP1_shRNA,
                                         hek_kd_mats,
                                         "HEK293 PTBP1 KD vs. Full Length Rescue", "HEK293 PTBP1 KD vs. E9 KO Rescue",
                                         slope = T)

# LAML analysis
laml_mats <- data.table::fread("~/transfer/SF_mutation/TCGA-LAML.SE.U2AF1_S34F_sorted.MATS.df.txt")

laml_sparks_result <- data.table::fread("/Users/harryyang/transfer/SF_Mutation/TCGA-LAML_result.min_0.txt")
# laml_sparks_result <- perform_SPARKS_analysis_with_overlap_filter(laml_mats,
#                                                                   kd_library,
#                                                                   "TCGA-LAML_U2AF1_S34F",
#                                                                   num_cores = 3)
#
# # add the S34F data
# laml_new_result <- add_custom_library_to_SPARKS_test_result(laml_mats,
#                                                             laml_sparks_result,
#                                                             k562_mats,
#                                                             "K562_U2AF1_S34F",
#                                                             kd_library)
# add_plot_title(generate_enrichment_barplot(laml_new_result,
#                                            bar_color = "violetred1",
#                                            num_plot = 10,
#                                            select_genes = c("U2AF1", "PTBP1")),
#                "TCGA-LAML U2AF1 S34F",
#                title_color = "violetred1") # 350 x 400
generate_enrichment_barplot(laml_sparks_result,
                            bar_color = "violetred1",
                            num_plot = 10,
                            select_genes = c("U2AF1", "PTBP1"),
                            manual_colors = list("PTBP1" = "darkorange2",
                                                 "U2AF1" = "violetred1")) +
  ggtitle("TCGA-LAML U2AF1 S34F") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

# generate heatmap with K562
generate_RBP_KD_correlation_scatter_plot(laml_mats,
                                         k562_mats,
                                         "TCGA-LAML U2AF1 S34F vs. WT",
                                         "Ilagan - K562 U2AF1 S34F vs. WT",
)

# generate heatmap with KD
k562_kd_sparks_file <- '/Users/harryyang/transfer/SF_Mutation/K562_U2AF1_KD.30.SPARKS.rds'
k562_kd_sparks <- readRDS(k562_kd_sparks_file)
k562_kd_mats <- import_SPARKS_MATS_for_analysis(k562_kd_sparks, "SE")

generate_RBP_KD_correlation_scatter_plot(laml_mats,
                                         k562_kd_mats,
                                         "TCGA-LAML U2AF1 S34F vs. WT",
                                         "Ilagan - K562 U2AF1 KD vs. WT",
)


generate_RBP_KD_correlation_scatter_plot(laml_mats,
                                         kd_library$K562_U2AF1_shRNA,
                                         "TCGA-LAML U2AF1 S34F vs. WT",
                                         "ENCODE - K562 U2AF1 KD vs. WT",
)


laml_hek_kd_result <- add_custom_library_to_SPARKS_test_result(laml_mats,
                                                               laml_sparks_result,
                                                               hek_kd_mats,
                                                               "HEK293_PTBP1_siRNA",
                                                               kd_library)
laml_hek_ko_result <- add_custom_library_to_SPARKS_test_result(laml_mats,
                                                               laml_hek_kd_result,
                                                               hek_e9koki_mats,
                                                               "HEK293_PTBP1_E9KO",
                                                               kd_library)

add_plot_title(generate_enrichment_barplot(laml_hek_ko_result,
                                           bar_color = "violetred1",
                                           num_plot = 10,
                                           select_genes = c("U2AF1", "PTBP1")),
               "TCGA-LAML U2AF1 S34F",
               title_color = "violetred1") # 350 x 400

##### BOXPLOTS #####

# general plot for example events
target_events <- c("FXR1:chr3:+:180688146:180693100:180693192:180693909:SE",
                   "PTBP1:chr19:+:805187:805491:805569:806407:SE",
                   "STRAP:chr12:+:16035753:16036474:16036610:16042861:SE",
                   "ATR:chr3:-:142168444:142169307:142169444:142171969:SE")

laml_target_entry <- subset(laml_mats, event %in% target_events)
laml_target_plot_df <- do.call(rbind, apply(laml_target_entry, 1, function(x) {
  # print(as.numeric(strsplit(x['psi_values'], ",")[[1]]))
  # print(do.call(as.numeric, strsplit(x['psi_values'], ",")))
  entry <- data.frame(psi = as.numeric(strsplit(x['psi_values'], ",")[[1]]),
                      genotype = c(rep("WT", length(laml_s34_wt_index)), rep("S34F", length(laml_s34_mut_index))),
                      gene = strsplit(x['event'], ":")[[1]][1])
  return(entry)
}))


# change order so WT could come first
laml_target_plot_df$genotype <- factor(laml_target_plot_df$genotype, levels = c("WT", "S34F"))


# generate plot
p_box_laml <- ggplot(laml_target_plot_df,
                     aes(x = genotype,
                         y = psi,
                         fill = genotype)) +
  geom_boxplot(outlier.shape = NA,
               width = 0.5) +
  geom_jitter(alpha = 0.2,
              width = 0.25) +
  scale_fill_manual(values = c("WT" = "royalblue3",
                               "S34F" = "lightskyblue1")) +
  labs(x = "U2AF1 Genotype in TCGA-LAML",
       y = "Percent Spliced In (PSI)") +
  scale_x_discrete(labels = c("WT" = sprintf("WT\n(n = %i)", length(laml_s34_wt_index)),
                              "S34F" = sprintf("S34F\n(n = %i)", length(laml_s34_mut_index)))) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1.2),
                     expand = expansion(mult = c(0, 0)),
                     breaks = seq(0, 4)/4) +
  # labs(subtitle = expression(p <= 3.7670*10^{-64})) +
  theme(axis.text.y   = element_text(size = 8),
        axis.text.x   = element_text(size = 8),
        axis.title.y  = element_text(size = 10),
        axis.title.x  = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)
  ) +
  ggsignif::geom_signif(comparisons = list(c("WT", "S34F")),
                        map_signif_level = TRUE,
                        size = 1,
                        margin_top = 0.1,
                        annotations = c("***"),
                        tip_length = 0.01) +
  facet_wrap(~gene, nrow = 1)

# calculate for K562

k562_mut_target_entry <- subset(k562_mats, event %in% target_events)
k562_mut_target_plot_df <- do.call(rbind, apply(k562_mut_target_entry, 1, function(x) {
  # print(as.numeric(strsplit(x['psi_values'], ",")[[1]]))
  # print(do.call(as.numeric, strsplit(x['psi_values'], ",")))
  entry <- data.frame(psi = as.numeric(strsplit(x['psi_values'], ",")[[1]]),
                      genotype = c(rep("WT", 2), rep("S34F", 1)),
                      gene = strsplit(x['event'], ":")[[1]][1])
  return(entry)
}))


# change order so WT could come first
k562_mut_target_plot_df$genotype <- factor(k562_mut_target_plot_df$genotype, levels = c("WT", "S34F"))


# generate plot
p_box_k562_mut <- ggplot(k562_mut_target_plot_df,
                         aes(x = genotype,
                             y = psi,
                             fill = genotype)) +
  geom_boxplot(outlier.shape = NA,
               width = 0.5) +
  geom_jitter(alpha = 0.2,
              width = 0.25) +
  scale_fill_manual(values = c("WT" = "royalblue3",
                               "S34F" = "lightskyblue1")) +
  labs(x = "U2AF1 Genotype in TCGA-LAML",
       y = "Percent Spliced In (PSI)") +
  # scale_x_discrete(labels = c("WT" = sprintf("WT\n(n = %i)", length(k562_mut_s34_wt_index)),
  # "S34F" = sprintf("S34F\n(n = %i)", length(k562_mut_s34_mut_index)))) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1.2),
                     expand = expansion(mult = c(0, 0)),
                     breaks = seq(0, 4)/4) +
  # labs(subtitle = expression(p <= 3.7670*10^{-64})) +
  theme(axis.text.y   = element_text(size = 8),
        axis.text.x   = element_text(size = 8),
        axis.title.y  = element_text(size = 10),
        axis.title.x  = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)
  ) +
  ggsignif::geom_signif(comparisons = list(c("WT", "S34F")),
                        map_signif_level = TRUE,
                        size = 1,
                        margin_top = 0.1,
                        annotations = c("***"),
                        tip_length = 0.01) +
  facet_wrap(~gene, nrow = 1)


# calculate for K562

k562_kd_target_entry <- subset(k562_kd_mats, event %in% target_events)
k562_kd_target_plot_df <- do.call(rbind, apply(k562_kd_target_entry, 1, function(x) {
  # print(as.numeric(strsplit(x['psi_values'], ",")[[1]]))
  # print(do.call(as.numeric, strsplit(x['psi_values'], ",")))
  entry <- data.frame(psi = as.numeric(strsplit(x['psi_values'], ",")[[1]]),
                      genotype = c(rep("WT", 2), rep("KD", 1)),
                      gene = strsplit(x['event'], ":")[[1]][1])
  return(entry)
}))


# change order so WT could come first
k562_kd_target_plot_df$genotype <- factor(k562_kd_target_plot_df$genotype, levels = c("WT", "KD"))


# generate plot
# p_box_k562_kd <-
  ggplot(k562_kd_target_plot_df,
                         aes(x = genotype,
                             y = psi,
                             fill = genotype)) +
  geom_boxplot(outlier.shape = NA,
               width = 0.5) +
  geom_jitter(alpha = 0.2,
              width = 0.25) +
  scale_fill_manual(values = c("WT" = "royalblue3",
                               "KD" = "lightskyblue1")) +
  labs(x = "U2AF1 Genotype in TCGA-LAML",
       y = "Percent Spliced In (PSI)") +
  # scale_x_discrete(labels = c("WT" = sprintf("WT\n(n = %i)", length(k562_kd_s34_wt_index)),
  # "S34F" = sprintf("S34F\n(n = %i)", length(k562_kd_s34_kd_index)))) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1.2),
                     expand = expansion(mult = c(0, 0)),
                     breaks = seq(0, 4)/4) +
  # labs(subtitle = expression(p <= 3.7670*10^{-64})) +
  theme(axis.text.y   = element_text(size = 8),
        axis.text.x   = element_text(size = 8),
        axis.title.y  = element_text(size = 10),
        axis.title.x  = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)
  ) +
  ggsignif::geom_signif(comparisons = list(c("WT", "KD")),
                        map_signif_level = TRUE,
                        size = 1,
                        margin_top = 0.1,
                        # test = 't.test',

                        annotations = c("***", ">", "?>",":"),
                        tip_length = 0.01) +
  facet_wrap(~gene, nrow = 1)



