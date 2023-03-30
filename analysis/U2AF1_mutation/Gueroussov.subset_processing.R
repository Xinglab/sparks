##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(fgsea)
library(metap)
library(dplyr)

input_kd_library <- readRDS('/home/yangt3/xinglab/SPARKS/library/Clean_minfiltered.SE.library.rds')

spl_types <- c("SE", "A3SS", "A5SS")
num_cores <- 8
overlap_percent <- 30
overlap_ratio <- overlap_percent / 100



input_sparks <- readRDS('/home/yangt3/xinglab/SPARKS/final_results/Guerrosov_PTBP12_siRNA.SPARKS.rds')

# generate subset
colnames(input_sparks@psi_df$SE)


# generate subset for others
kd_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                          input_kd_library,
                                          "HEK293_PTBP1_KD",
                                          c(5, 6),
                                          c(1, 2),
                                          num_cores = num_cores,
                                          overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(kd_subset, sprintf("HEK293_PTBP1_KD.%s.SPARKS.rds", overlap_percent))


# generate subset for others
fl_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                          input_kd_library,
                                          "HEK293_PTBP1_FLrescue",
                                          c(7, 8),
                                          c(1, 2),
                                          num_cores = num_cores,
                                          overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(fl_subset, sprintf("HEK293_PTBP1_FLrescue.%s.SPARKS.rds", overlap_percent))

# generate subset for others
e9_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                          input_kd_library,
                                          "HEK293_PTBP1_delE9rescue",
                                          c(3, 4),
                                          c(1, 2),
                                          num_cores = num_cores,
                                          overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(e9_subset, sprintf("HEK293_PTBP1_delE9rescue.%s.SPARKS.rds", overlap_percent))

# generate subset for others
dele9_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                             input_kd_library,
                                             "HEK293_PTBP1_delE9",
                                             c(5, 6),
                                             c(3, 4),
                                             num_cores = num_cores,
                                             overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(dele9_subset, sprintf("HEK293_PTBP1_delE9.%s.SPARKS.rds", overlap_percent))
