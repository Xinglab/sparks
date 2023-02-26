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



input_sparks <- readRDS('./Fei_LUAD_U2AF1_S34F.SPARKS.rds')

# generate subset
colnames(input_sparks@psi_df$SE)

h441_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                            input_kd_library,
                                            "H441_U2AF1_S34F",
                                            c(1, 2),
                                            c(7, 8),
                                            num_cores = num_cores,
                                            overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(h441_subset, sprintf("H441_U2AF1_S34F.%s.SPARKS.rds", overlap_percent))


hbec3kt_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                            input_kd_library,
                                            "HBEC3kt_U2AF1_S34F",
                                            c(3, 4),
                                            c(9, 10),
                                            num_cores = num_cores,
                                            overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(hbec3kt_subset, sprintf("HBEC3kt_U2AF1_S34F.%s.SPARKS.rds", overlap_percent))


hcc78_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                            input_kd_library,
                                            "HCC78_U2AF1_S34F",
                                            c(5, 6),
                                            c(11, 12),
                                            num_cores = num_cores,
                                            overlap_ratio_threshold = overlap_ratio)

# save the subset sparks
saveRDS(hcc78_subset, sprintf("HCC78_U2AF1_S34F.%s.SPARKS.rds", overlap_percent))
