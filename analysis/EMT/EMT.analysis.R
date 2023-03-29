##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(dplyr)

##### FUNCTIONS #####



kd_library <- readRDS("/Users/harryyang/transfer/Clean_minfiltered.SE.library.rds")

##### ESRP ANALYSIS #####
## test similarity between ESRP1/2 KD experiemnts
# read in LNCaP ESRP1/2 siRNA data
lncap_esrp_sparks_file <- '/Users/harryyang/transfer/EMT/Munkley_LNCaP_ESRP12_siRNA.SPARKS.rds'
lncap_esrp_sparks <- readRDS(lncap_esrp_sparks_file)
lncap_esrp_mats <- import_SPARKS_MATS_for_analysis(lncap_esrp_sparks, "SE")

lncap_esrp_result <- lncap_esrp_sparks@SPARKS_analysis_result$SE

# read in H358 ESRP1/2 siRNA data
h358_esrp_sparks_file <- '/Users/harryyang/transfer/EMT/Yang_ESRP_KD.SPARKS.rds'
h358_esrp_sparks <- readRDS(h358_esrp_sparks_file)
h358_esrp_mats <- import_SPARKS_MATS_for_analysis(h358_esrp_sparks, "SE")

h358_esrp_result <- h358_esrp_sparks@SPARKS_analysis_result$SE

# calculate and add the cross test result
lncap_new_result <- add_custom_library_to_SPARKS_test_result(lncap_esrp_mats,
                                                             lncap_esrp_result,
                                                             h358_esrp_mats,
                                                             "H358_ESRP1/2_siRNA",
                                                             kd_library)
h358_new_result <- add_custom_library_to_SPARKS_test_result(h358_esrp_mats,
                                                            h358_esrp_result,
                                                            lncap_esrp_mats,
                                                            "LNCaP_ESRP1/2_siRNA",
                                                            kd_library)


# generate scatter plot
generate_RBP_KD_correlation_scatter_plot(lncap_esrp_mats,
                                         h358_esrp_mats,
                                         "LNCaP ESRP1/2 siRNA vs. Control",
                                         "H358 ESRP1/2 siRNA vs. Control") +
  theme(axis.title  = element_text(color = "dodgerblue2"))

# generate bar plot
generate_enrichment_barplot(lncap_new_result,
                            bar_color = "dodgerblue1",
                            num_plot = 10,
                            manual_colors = list("ESRP1/2" = "dodgerblue3")) +
  ggtitle("LNCaP ESRP1/2 siRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold",
                                  color = "dodgerblue2"))
generate_enrichment_barplot(h358_new_result,
                            bar_color = "dodgerblue1",
                            num_plot = 10,
                            manual_colors = list("ESRP1/2" = "dodgerblue3")) +
  ggtitle("H358 ESRP1/2 siRNA") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold",
                                  color = "dodgerblue2"))


##### ZEB1 ANALYSIS #####
### Yang H358 analysis
# read in the data
yang_sparks <- readRDS('/Users/harryyang/transfer/EMT/Yang_H358_EMT.SPARKS.rds')
yang_mats <- import_SPARKS_MATS_for_analysis(yang_sparks, "SE")
yang_result_basic <- yang_sparks@SPARKS_analysis_result$SE

# perform custom SPARKS analysis
yang_result_h358 <- add_custom_library_to_SPARKS_test_result(yang_mats,
                                                             yang_result_basic,
                                                             h358_esrp_mats,
                                                             "H358_ESRP1/2_siRNA",
                                                             kd_library)
yang_result_esrp <- add_custom_library_to_SPARKS_test_result(yang_mats,
                                                             yang_result_h358,
                                                             lncap_esrp_mats,
                                                             "LNCaP_ESRP1/2_siRNA",
                                                             kd_library)

# generate enrichment plot
generate_enrichment_barplot(yang_result_esrp,
                            bar_color = "indianred2",
                            num_plot = 10,
                            manual_colors = list("ESRP1/2" = "dodgerblue3")) +
  ggtitle("H358 ZEB1 OE") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold",
                                  color = "indianred2"))



### Zhang HCC827 analysis
# perform SPARKS analysis
zhang_sparks <- readRDS('/Users/harryyang/transfer/EMT/Zhang_HCC827_ZEB1_OE.SPARKS.rds')
zhang_mats <- import_SPARKS_MATS_for_analysis(zhang_sparks, "SE")
zhang_result_basic <- zhang_sparks@SPARKS_analysis_result$SE

zhang_result_h358 <- add_custom_library_to_SPARKS_test_result(zhang_mats,
                                                              zhang_result_basic,
                                                              h358_esrp_mats,
                                                              "H358_ESRP1/2_siRNA",
                                                              kd_library)
zhang_result_esrp <- add_custom_library_to_SPARKS_test_result(zhang_mats,
                                                              zhang_result_h358,
                                                              lncap_esrp_mats,
                                                              "LNCaP_ESRP1/2_siRNA",
                                                              kd_library)

generate_enrichment_barplot(zhang_result_esrp,
                                           bar_color = "indianred2",
                                           num_plot = 10,
                            manual_colors = list("ESRP1/2" = "dodgerblue3")) +
  ggtitle("HCC827 ZEB1 OE") +
  theme(plot.title = element_text(size = 10,
                                  hjust = 0.5,
                                  face = "bold",
                                  color = "indianred2"))




##### EXPRESSION #####

# generate expression plot for ZEB1 and ESRP1
yang_exp_df <- query_expression_data_for_gene_set(yang_sparks, c("ESRP1", "ESRP2", "ZEB1"))

yang_exp_melt <- reshape2::melt(yang_exp_df)

# add relevant information
yang_exp_melt$status <- unlist(lapply(yang_exp_melt$variable, function(x) strsplit(as.character(x), "_")[[1]][2]))
yang_exp_melt$condition <- ifelse(yang_exp_melt$status == "No",  # No means no ZEB1 Over expression
                                  "Control",
                                  "Overexpression")
p_yang_exp <- ggplot(yang_exp_melt,
                     aes(x = condition,
                         y = log10(value + 1))) +
  geom_boxplot() +
  geom_jitter(aes(color = condition)) +
  facet_grid(~ gene) +
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
  )+
  ggsignif::geom_signif(comparisons = list(c("Control", "Overexpression")),
                        # map_signif_level = TRUE,
                        size = 1,
                        margin_top = 0.1,
                        test = "t.test",
                        tip_length = 0.01,
                        annotations = c("***")) +
  labs(x = "H358 ZEB1 Status",
       y = "Expression in log10(TPM + 1)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "paleturquoise3",
                                "Overexpression" = "royalblue")) +
  scale_y_continuous(limits = c(0, 4),
                     expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5))





# generate expression plot for ZEB1 and ESRP1
zhang_exp_df <- query_expression_data_for_gene_set(zhang_sparks, c("ESRP1", "ESRP2", "ZEB1"))

zhang_exp_melt <- reshape2::melt(zhang_exp_df)

# add relevant information
zhang_exp_melt$status <- unlist(lapply(zhang_exp_melt$variable, function(x) strsplit(as.character(x), "_")[[1]][2]))
zhang_exp_melt$condition <- ifelse(zhang_exp_melt$status == "Control",  # No means no ZEB1 Over expression
                                   "Control",
                                   "Overexpression")
p_zhang_exp <- ggplot(zhang_exp_melt,
                      aes(x = condition,
                          y = log10(value + 1))) +
  geom_boxplot() +
  geom_jitter(aes(color = condition)) +
  facet_grid(~ gene) +
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
  )+
  ggsignif::geom_signif(comparisons = list(c("Control", "Overexpression")),
                        size = 1,
                        margin_top = 0.48,
                        test = "t.test",
                        tip_length = 0.01,
                        annotations = c("***")) +
  labs(x = "HCC827 ZEB1 Status",
       y = "Expression in log10(TPM + 1)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "paleturquoise3",
                                "Overexpression" = "royalblue")) +
  scale_y_continuous(limits = c(0, 4),
                     expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5))

cowplot::plot_grid(p_yang_exp, p_zhang_exp, align = "h", axis = "tb")




