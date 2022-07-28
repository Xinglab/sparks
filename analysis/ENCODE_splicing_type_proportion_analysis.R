##### CODE #####
library(SPARKS)
library(ggplot2)
library(dplyr)


##### FUNCTIONS #####
generate_splicing_summary_table <- function(study_mats_list, study){

  test_sig_df  <- do.call(rbind, lapply(names(study_mats_list), function(spl_type){
    test_mats <- study_mats_list[[spl_type]]
    test_sig_subset <- subset(test_mats, abs(beta) > 0.1)
    test_sig_subset$spl_type <- spl_type

    return(test_sig_subset)
  }))

  ## Directional
  # test_sig_count_df  <- data.frame(table(data.frame(ifelse(test_sig_df $beta > 0, "pos", "neg"), test_sig_df $spl_type)))
  # colnames(test_sig_count_df ) <- c("direction", "spl_type", "freq")
  # test_sig_count_df$loc <- ifelse(test_sig_count_df $direction == "neg", max(test_sig_df $beta) + 0.1, min(test_sig_df $beta) - 0.1)

  ## Non-directional
  test_sig_count_df <- data.frame(table(data.frame(abs(test_sig_df$beta) >= 0, test_sig_df $spl_type)))
  colnames(test_sig_count_df) <- c("direction", "spl_type", "freq")
  test_sig_count_df$direction <- NULL
  test_sig_count_df$study <- study
  return(test_sig_count_df)
}

prepare_library_data_for_splicing_sumamry <- function(kd_library_all, study){

  output_list <- list()
  spl_types <- c("SE", "A3SS", "A5SS")

  dummy <- lapply(spl_types, function(spl_type){
    study_mats <- kd_library_all[[spl_type]][[study]]
    output_list[[spl_type]] <<- study_mats
  })
  return(output_list)
}


# import library
merged_kd_library <- readRDS('/Users/harryyang/research/Xing_Lab/sparks/data/ENCODE_library.KD_and_KO.count_20.rds')



# extract library names for combined summary
exp_names <- names(merged_kd_library$SE)

# run thru and get summary
splicing_summary_data <- do.call(rbind, lapply(exp_names, function(experiment){
  print(experiment)

  spl_data <- prepare_library_data_for_splicing_sumamry(merged_kd_library, experiment)
  summary_data <- generate_splicing_summary_table(spl_data, experiment)
  return(summary_data)
}))

# get summary for sorting
total_num_event_df <- as.data.frame(splicing_summary_data %>% group_by(study) %>% summarize(sum = sum(freq)) %>% arrange(-sum))
total_num_event_df$study <- factor(total_num_event_df$study, levels = total_num_event_df$study)

splicing_summary_data$study <- factor(splicing_summary_data$study, levels = total_num_event_df$study)

# generate plot
p_prop <- ggplot(splicing_summary_data, aes(x = study, y = freq, fill = spl_type)) +
  geom_bar(stat = "identity", position = "fill", width = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = "bottom") +
  labs(x = "ENCODE RBP Perturbation Experiments\n(sorted by total number of events)",
       y = "Proportion of AS events",
       fill = "Splicing Type") +
  scale_fill_manual(values = c("SE" = "darkgoldenrod1",
                               "A3SS" = "peachpuff3",
                               "A5SS" = "palevioletred1")) +
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
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0, 0))

p_num <- ggplot(total_num_event_df, aes(x = study, y = sum)) +
  geom_bar(stat = "identity", width = 1, fill = "olivedrab")+
  labs(x = "ENCODE RBP Perturbation Experiments\n(sorted by total number of events)",
       y = "Total Number of AS events")+
  # theme_minima
  theme(axis.text.y   = element_text(size=8),
        axis.text.x   = element_text(size=8),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.position = "bottom")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0, 0))

cowplot::plot_grid(p_num +
                     theme(axis.title.x = element_blank()),
                   p_prop + theme(legend.position = "none"),
                   cowplot::get_legend(p_prop),
                   align = 'hv',
                   axis = "tblr",
                   ncol = 1,
                   rel_heights = c(4, 4, 1)) ## size - w350 x h700
