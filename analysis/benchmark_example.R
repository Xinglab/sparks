##### CODE #####
library(START)
library(ggplot2)
score_method <- "GSEA"
spl_types <- c("SE", "A3SS", "A5SS")



##### ANALYSIS #####
merged_kd_library <- readRDS('/Users/harryyang/research/Xing_Lab/START/analysis/library/ENCODE.KDKO.signatures.MATS_df_list.count_20.rds')

# determine replicate RBP experiments
input_study_list <- names(merged_kd_library$SE)[grep("shRNA", names(merged_kd_library$SE))]
# generate study -> RBP table
rbp_info_df <- data.frame(study = input_study_list,
                          rbp = unlist(lapply(input_study_list, function(x) strsplit(x, "_")[[1]][2])))
rbp_num_exp <- table(rbp_info_df$rbp)
rbp_kd2_list <- names(rbp_num_exp)[rbp_num_exp == 2]

encode_rep_experiments <- rbp_info_df[rbp_info_df$rbp %in% rbp_kd2_list, ]$study

test_result_list <- lapply(spl_types, function(spl_type){
  input_mats <- merged_kd_library[[spl_type]][["K562_SRSF1_shRNA"]]

  hepg2_library <- merged_kd_library[[spl_type]][encode_rep_experiments[grep("HepG2", encode_rep_experiments)]]

  test_result <- perform_SPARKS_analysis(input_mats, hepg2_library, "K562_SRSF1_shRNA")

  return(test_result)

})
names(test_result_list) <- spl_types

add_plot_title(generate_strip_lollipop_plot_vertical(test_result_list, "SRSF1", "SE",
                                                     y_axis_title = "HepG2 Enrichment Score Rank",
                                                     strip_break_interval = 50),
               "K562 SRSF1 shRNA")

