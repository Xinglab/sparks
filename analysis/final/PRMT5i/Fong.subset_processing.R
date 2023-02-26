##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(fgsea)
library(metap)
library(dplyr)

#### FUNCTIONS #####
input_kd_library <- readRDS('/home/yangt3/xinglab/SPARKS/library/Clean_minfiltered.SE.library.rds')

spl_types <- c("SE", "A3SS", "A5SS")
num_cores <- 8
overlap_percent <- 30
overlap_ratio <- overlap_percent / 100



input_sparks <- readRDS('/home/yangt3/xinglab/SPARKS/final_results/PRMT5i/Fong_PRMT5i_2019.SPARKS.rds')

# generate subset
colnames(input_sparks@psi_df$SE)

# # add the S34F library
# snrpb_sparks <- readRDS('/home/yangt3/xinglab/START/snrpb/Correa_et_al_2016/Correa_SNRPB_siRNA.SPARKS.rds')
#
# snrpb_mats <- import_SPARKS_MATS_for_rerun(snrpb_sparks)
#
#
# # change the event
# input_kd_library[["U251_SNRPB_siRNA"]] <- snrpb_mats


wt_subset <- generate_subset_SPARKS_rerun(input_sparks,
                                            input_kd_library,
                                            "K562_PRMT5_inhibition",
                                            c(1, 2, 3),
                                            c(7, 8, 9),
                                            num_cores = num_cores,
                                            overlap_ratio_threshold = overlap_ratio)


# save the subset sparks
saveRDS(wt_subset, sprintf("K562_PRMT5_inhibition.%s.SPARKS.rds", overlap_percent))


