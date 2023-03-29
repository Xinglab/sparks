srsf1_original_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.MCF7_SRSF1_shRNA.SE.30.original.txt')
p_srsf1 <- generate_enrichment_barplot(srsf1_original_result_df,
                                       num_plot = 10,
                                       bar_color = "indianred1",
                                       manual_colors = list("SRSF1" = "red")) +
  ggtitle("MCF7 SRSF1 shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))





ptbp1_original_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.HEK293_PTBP1_KD.SE.30.original.txt')
p_ptbp1 <- generate_enrichment_barplot(ptbp1_original_result_df,
                                       num_plot = 10,
                                       bar_color = "darkorange",
                                       manual_colors = list("PTBP1" = "red")) +
  ggtitle("HEK293 PTBP1/2 siRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

srsf3_original_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.MDA-MB231_SRSF3_shRNA.SE.30.original.txt')
p_srsf3 <- generate_enrichment_barplot(srsf3_original_result_df,
                                       num_plot = 10,
                                       bar_color = "steelblue1",
                                       manual_colors = list("SRSF3" = "red")) +
  ggtitle("MDA-MB-231 SRSF3 shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

hnrnpk_original_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.Epiderm_HNRNPK_shRNA.SE.30.original.txt')
p_hnrnpk <- generate_enrichment_barplot(hnrnpk_original_result_df,
                                        num_plot = 10,
                                        bar_color = "thistle3",
                                        manual_colors = list("HNRNPK" = "red")) +
  ggtitle("Epiderm HNRNPK shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

cowplot::plot_grid(p_srsf1,
                   p_ptbp1,
                   p_srsf3,
                   p_hnrnpk,
                   nrow = 2) # 700 x 800

##### CELL LINE SWAP FIGURES #####
# SRSF1
srsf1_hepg2_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.MCF7_SRSF1_shRNA.SE.30.HepG2_swapped.txt')
p_srsf1_hepg2 <- generate_enrichment_barplot(srsf1_hepg2_result_df,
                                             num_plot = 10,
                                             bar_color = "indianred1",
                                             manual_colors = list("SRSF1" = "red")) +
  ggtitle("HepG2 SRSF1 shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

srsf1_k562_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.MCF7_SRSF1_shRNA.SE.30.K562_swapped.txt')
p_srsf1_k562 <- generate_enrichment_barplot(srsf1_k562_result_df,
                                            num_plot = 10,
                                            bar_color = "indianred1",
                                            manual_colors = list("SRSF1" = "red")) +
  ggtitle("K562 SRSF1 shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

cowplot::plot_grid(p_srsf1,
                   p_ptbp1,
                   p_srsf3,
                   p_hnrnpk,
                   p_srsf1_hepg2,
                   p_srsf1_k562,
                   nrow = 3) # 700 x 1200

# PTBP1
ptbp1_hepg2_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.HEK293_PTBP1_KD.SE.30.HepG2_swapped.txt')
p_ptbp1_hepg2 <- generate_enrichment_barplot(ptbp1_hepg2_result_df,
                                             num_plot = 10,
                                             bar_color = "darkorange",
                                             manual_colors = list("PTBP1" = "red")) +
  ggtitle("HepG2 PTBP1 shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

ptbp1_k562_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.HEK293_PTBP1_KD.SE.30.K562_swapped.txt')
p_ptbp1_k562 <- generate_enrichment_barplot(ptbp1_k562_result_df,
                                            num_plot = 10,
                                            bar_color = "darkorange",
                                            manual_colors = list("PTBP1" = "red")) +
  ggtitle("K562 PTBP1 shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

# SRSF3
srsf3_hepg2_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.MDA-MB231_SRSF3_shRNA.SE.30.HepG2_swapped.txt')
p_srsf3_hepg2 <- generate_enrichment_barplot(srsf3_hepg2_result_df,
                                             num_plot = 10,
                                             bar_color = "steelblue1",
                                             manual_colors = list("SRSF3" = "red")) +
  ggtitle("HepG2 SRSF3 shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

srsf3_k562_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.MDA-MB231_SRSF3_shRNA.SE.30.K562_swapped.txt')
p_srsf3_k562 <- generate_enrichment_barplot(srsf3_k562_result_df,
                                            num_plot = 10,
                                            bar_color = "steelblue1",
                                            manual_colors = list("SRSF3" = "red")) +
  ggtitle("K562 SRSF3 CRISPR") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

# HNRNPK
hnrnpk_hepg2_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.Epiderm_HNRNPK_shRNA.SE.30.HepG2_swapped.txt')
p_hnrnpk_hepg2 <- generate_enrichment_barplot(hnrnpk_hepg2_result_df,
                                              num_plot = 10,
                                              bar_color = "thistle3",
                                              manual_colors = list("HNRNPK" = "red")) +
  ggtitle("HepG2 HNRNPK shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))

hnrnpk_k562_result_df <- data.table::fread('/Users/harryyang/transfer/single_perturbation/benchmark_result.Epiderm_HNRNPK_shRNA.SE.30.K562_swapped.txt')
p_hnrnpk_k562 <- generate_enrichment_barplot(hnrnpk_k562_result_df,
                                             num_plot = 10,
                                             bar_color = "thistle3",
                                             manual_colors = list("HNRNPK" = "red")) +
  ggtitle("K562 HNRNPK shRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold"))


# combine the plots
cowplot::plot_grid(p_ptbp1_hepg2, p_srsf3_hepg2, p_hnrnpk_hepg2,
                   p_ptbp1_k562, p_srsf3_k562, p_hnrnpk_k562,
                   nrow = 2) # 1050 x 800
