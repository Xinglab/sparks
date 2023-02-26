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



input_sparks <- readRDS('./Liu_BrCa_SF3B1_K700E.SPARKS.rds')

# generate subset
colnames(input_sparks@psi_df$SE)

mcf10_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                            input_kd_library,
                                            "MCF10A_SF3B1_K700E",
                                            c(1, 2),
                                            c(7, 8),
                                            num_cores = num_cores,
                                            overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(mcf10_subset, sprintf("MCF10A_SF3B1_K700E.%s.SPARKS.rds", overlap_percent))


mcf7_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                               input_kd_library,
                                               "MCF7_SF3B1_K700E",
                                               c(3, 4),
                                               c(9, 10),
                                               num_cores = num_cores,
                                               overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(mcf7_subset, sprintf("MCF7_SF3B1_K700E.%s.SPARKS.rds", overlap_percent))


t47d_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                             input_kd_library,
                                             "T47D_SF3B1_K700E",
                                             c(5, 6),
                                             c(11, 12),
                                             num_cores = num_cores,
                                             overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(t47d_subset, sprintf("T47D_SF3B1_K700E.%s.SPARKS.rds", overlap_percent))
