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



input_sparks <- readRDS('./Ilagan_U2AF1_mutant.SPARKS.rds')

# generate subset
colnames(input_sparks@psi_df$SE)

s34f_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                              input_kd_library,
                                              "K562_U2AF1_S34F",
                                              c(1, 2),
                                              c(5),
                                              num_cores = num_cores,
                                              overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(s34f_subset, sprintf("K562_U2AF1_S34F.%s.SPARKS.rds", overlap_percent))

# generate subset for others
kd_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                            input_kd_library,
                                            "K562_U2AF1_KD",
                                            c(1, 2),
                                            c(7),
                                            num_cores = num_cores,
                                            overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(kd_subset, sprintf("K562_U2AF1_KD.%s.SPARKS.rds", overlap_percent))

