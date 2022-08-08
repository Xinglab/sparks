##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(dplyr)

##### ANALYSIS #####
# import library
merged_kd_library <- readRDS('/Users/harryyang/Documents/Research/Xing_Lab/sparks/data/ENCODE_library.KD_and_KO.count_20.rds')
spl_types <- c("SE", "A3SS", "A5SS")
score_method <- "GSEA"

##### ESRP ANALYSIS #####
## test similarity between ESRP1/2 KD experiemnts
# read in LNCaP ESRP1/2 siRNA data
lncap_esrp_sparks_file <- '/Users/harryyang/research/Xing_Lab/sparks/data/EMT/Munkley_LNCaP_ESRP12_siRNA.SPARKS.rds'
lncap_esrp_start <- readRDS(lncap_esrp_sparks_file)
lncap_esrp_study <- "LNCaP_ESRP1/2_siRNA"
lncap_esrp_mats <- import_SPARKS_MATS_for_analysis(lncap_esrp_start, "SE")

# read in H358 ESRP1/2 siRNA data
h358_esrp_sparks_file <- '/Users/harryyang/research/Xing_Lab/sparks/data/EMT/Yang_ESRP_KD.SPARKS.rds'
h358_esrp_start <- readRDS(h358_esrp_sparks_file)
h358_esrp_study <- "H358_ESRP1/2_siRNA"
h358_esrp_mats <- import_SPARKS_MATS_for_analysis(h358_esrp_start, "SE")

# generate scatter plot
generate_RBP_KD_correlation_scatter_plot(lncap_esrp_mats,
                                         h358_esrp_mats,
                                         "LNCaP ESRP1/2 siRNA vs. Control",
                                         "H358 ESRP1/2 siRNA vs. Control") +
  theme(axis.title  = element_text(color = "dodgerblue2"))

# import H358 data as library for LNCaP analysis
h358_esrp_library <- import_custom_study_as_SPARKS_library(h358_esrp_sparks_file,
                                                           h358_esrp_study)
# import LNCaP data as library for H358 analysis
lncap_esrp_library <- import_custom_study_as_SPARKS_library(lncap_esrp_sparks_file,
                                                            lncap_esrp_study)

# run SPARKS using the custom ESRP library
lncap_esrp_result_list <- perform_SPARKS_analysis_for_all_splice_types(lncap_esrp_start,
                                                                       h358_esrp_library,
                                                                       test_study = lncap_esrp_study)
h358_esrp_result_list <- perform_SPARKS_analysis_for_all_splice_types(h358_esrp_start,
                                                                      lncap_esrp_library,
                                                                      test_study = h358_esrp_study)

# generate merged custom ESRP library for downstream run
esrp_library <- merge_custom_SPARKS_libraries(list(h358_esrp_library,
                                                   lncap_esrp_library))

# add the data to regular SPARKS result inside the object
h358_start_list <- generate_SPARKS_result_with_custom_library_results(h358_esrp_start,
                                                                      h358_esrp_result_list)
lncap_start_list <- generate_SPARKS_result_with_custom_library_results(lncap_esrp_start,
                                                                       lncap_esrp_result_list)

# generate lollipop plot for ESRP enrichment - w250 x h200
add_plot_title(generate_strip_lollipop_plot_vertical(lncap_start_list,
                                                     c("ESRP1/2"),
                                                     "SE",
                                                     manual_colors = "dodgerblue2"),
               "LNCaP ESRP1/2 siRNA",
               title_color = "dodgerblue2")
add_plot_title(generate_strip_lollipop_plot_vertical(h358_start_list,
                                                     c("ESRP1/2"),
                                                     "SE",
                                                     manual_colors = "dodgerblue2"),
               "H358 ESRP1/2 siRNA",
               title_color = "dodgerblue2")

# make supplementary plots - w250 x h400
add_plot_title(cowplot::plot_grid(generate_strip_lollipop_plot_vertical(lncap_start_list,
                                                                        c("ESRP1/2"),
                                                                        "A5SS",
                                                                        manual_colors = "dodgerblue2"),
                                  generate_strip_lollipop_plot_vertical(lncap_start_list,
                                                                        c("ESRP1/2"),
                                                                        "A3SS",
                                                                        manual_colors = "dodgerblue2"),
                                  ncol = 1),
               "LNCaP ESRP1/2 siRNA",
               title_color = "dodgerblue2")
add_plot_title(cowplot::plot_grid(generate_strip_lollipop_plot_vertical(h358_start_list,
                                                                        c("ESRP1/2"),
                                                                        "A5SS",
                                                                        manual_colors = "dodgerblue2"),
                                  generate_strip_lollipop_plot_vertical(h358_start_list,
                                                                        c("ESRP1/2"),
                                                                        "A3SS",
                                                                        manual_colors = "dodgerblue2"),
                                  ncol = 1),
               "H358 ESRP1/2 siRNA",
               title_color = "dodgerblue2")



##### EMT ANALYSIS #####
# perform SPARKS analysis
yang_start <- readRDS('/Users/harryyang/research/Xing_Lab/sparks/data/EMT/Yang_H358_EMT.SPARKS.rds')
yang_study <- "H358_ZEB1_OE"

yang_esrp_list <- perform_SPARKS_analysis_for_all_splice_types(yang_start,
                                                               esrp_library,
                                                               yang_study)
yang_result_list <- generate_SPARKS_result_with_custom_library_results(yang_start,
                                                                       yang_esrp_list)

add_plot_title(generate_strip_lollipop_plot_vertical(yang_result_list,
                                                     c("ESRP1/2"),
                                                     "SE",
                                                     manual_color = "dodgerblue2"),
               title_color = "indianred2",
               plot_title = "H358 ZEB1 OE")

add_plot_title(cowplot::plot_grid(generate_strip_lollipop_plot_vertical(yang_result_list,
                                                                        c("ESRP1/2"),
                                                                        "A5SS",
                                                                        manual_colors = "dodgerblue2"),
                                  generate_strip_lollipop_plot_vertical(yang_result_list,
                                                                        c("ESRP1/2"),
                                                                        "A3SS",
                                                                        manual_colors = "dodgerblue2"),
                                  ncol = 1),
               "H358 ZEB1 OE",
               title_color = "indianred2")

# generate expression plot for ZEB1 and ESRP1
yang_exp_df <- query_expression_data_for_gene_set(yang_start, c("ESRP1", "ESRP2", "ZEB1"))

yang_exp_melt <- reshape2::melt(yang_exp_df)

# add relevant information
yang_exp_melt$status <- unlist(lapply(yang_exp_melt$variable, function(x) strsplit(as.character(x), "_")[[1]][2]))
yang_exp_melt$condition <- ifelse(yang_exp_melt$status == "No",  # No means no ZEB1 Over expression
                                  "Control",
                                  "Overexpression")
ggplot(yang_exp_melt,
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
                        map_signif_level = TRUE,
                        size = 1,
                        margin_top = 0.1,
                        test = "t.test",
                        tip_length = 0.01) +
  labs(x = "H358 ZEB1 Status",
       y = "Expression in log10(TPM + 1)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "paleturquoise3",
                               "Overexpression" = "royalblue")) +
  scale_y_continuous(limits = c(0, 4),
                     expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5))




##### Other EMT data set #####
zhang_start <- readRDS('/Users/harryyang/research/Xing_Lab/sparks/data/EMT/Zhang_HCC827_ZEB1_OE.SPARKS.rds')
zhang_study <- "HCC827_ZEB1_OE"


zhang_esrp_list <- perform_SPARKS_analysis_for_all_splice_types(zhang_start,
                                                                esrp_library,
                                                                zhang_study)
zhang_result_list <- generate_SPARKS_result_with_custom_library_results(zhang_start,
                                                                        zhang_esrp_list)

add_plot_title(generate_strip_lollipop_plot_vertical(zhang_result_list,
                                                     c("ESRP1/2"),
                                                     "SE",
                                                     manual_color = "dodgerblue2"),
               title_color = "indianred2",
               plot_title = "HCC827 ZEB1 OE")

add_plot_title(cowplot::plot_grid(generate_strip_lollipop_plot_vertical(zhang_result_list, c("ESRP1/2"), "A5SS", manual_colors = "dodgerblue2"),
                                  generate_strip_lollipop_plot_vertical(zhang_result_list, c("ESRP1/2"), "A3SS", manual_colors = "dodgerblue2"),
                                  ncol = 1),
               "HCC827 ZEB1 OE", title_color = "indianred2")


# generate expression plot for ZEB1 and ESRP1
zhang_exp_df <- query_expression_data_for_gene_set(zhang_start, c("ESRP1", "ESRP2", "ZEB1"))

zhang_exp_melt <- reshape2::melt(zhang_exp_df)

# add relevant information
zhang_exp_melt$status <- unlist(lapply(zhang_exp_melt$variable, function(x) strsplit(as.character(x), "_")[[1]][2]))
zhang_exp_melt$condition <- ifelse(zhang_exp_melt$status == "Control",  # No means no ZEB1 Over expression
                                   "Control",
                                   "Overexpression")
ggplot(zhang_exp_melt,
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
                        map_signif_level = TRUE,
                        size = 1,
                        margin_top = 0.1,
                        test = "t.test",
                        tip_length = 0.01) +
  labs(x = "HCC827 ZEB1 Status",
       y = "Expression in log10(TPM + 1)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Control" = "paleturquoise3",
                                "Overexpression" = "royalblue")) +
  scale_y_continuous(limits = c(0, 4),
                     expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5))


##### Other EMT dataset #####
balest_start <- readRDS('/Users/harryyang/research/Xing_Lab/sparks/data/EMT/Balestrieri_ZEB1_KO.SPARKS.rds')
balest_study <- "MiaPaCa_ZEB1_KO"

balest_esrp_list <- perform_SPARKS_analysis_for_all_splice_types(balest_start,
                                                                 esrp_library,
                                                                 balest_study)
balest_result_list <- generate_SPARKS_result_with_custom_library_results(balest_start,
                                                                         balest_esrp_list)

add_plot_title(generate_strip_lollipop_plot_vertical(balest_result_list,
                                                     c("ESRP1/2"),
                                                     "SE",
                                                     manual_color = "dodgerblue2"),
               title_color = "indianred2",
               plot_title = "MiaPaCa2 ZEB1 KO")


add_plot_title(cowplot::plot_grid(generate_strip_lollipop_plot_vertical(balest_result_list,
                                                                        c("ESRP1/2"),
                                                                        "A5SS",
                                                                        manual_colors = "dodgerblue2"),
                                  generate_strip_lollipop_plot_vertical(balest_result_list,
                                                                        c("ESRP1/2"),
                                                                        "A3SS",
                                                                        manual_colors = "dodgerblue2"),
                                  ncol = 1),
               "MiaPaCa2 ZEB1 KO",
               title_color = "indianred2")

## Expression is not included because this is a single-ended read dataset
