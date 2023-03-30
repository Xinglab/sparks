##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(fgsea)
library(metap)
library(dplyr)


##### ANALYSIS CODE #####
# import library for analysis
print("Loading Database")
kd_library_all <- readRDS('/home/yangt3/xinglab/SPARKS/library/Clean_minfiltered.SE.library.rds')

# set up the parameters
output_dir <- "./"
count_threshold <- 20

spl_type = "SE"

# input the argument
args = commandArgs(trailingOnly=TRUE)
input_file <- args[1]
num_events <- as.integer(args[2])
num_ratio <- num_events / 100

# read in the file
input_spark <- readRDS(input_file)

input_study <- strsplit(strsplit(input_file, "/")[[1]][length(strsplit(input_file, "/")[[1]])], "\\.")[[1]][1]
print(input_study)

input_rbp <- strsplit(input_study, "_")[[1]][2]


# import study mats
study_mats <- import_SPARKS_MATS_for_rerun(input_spark, spl_type)


## remove each library for swap and add new
# find the library with the RBP KD - need to do this since SRSF3 is CRISPR actually
rbp_kd_lib_list <- names(kd_library_all)[grep(input_rbp, names(kd_library_all))]

# make new library without the specific experiment
hepg2_minus_library <- kd_library_all[names(kd_library_all) != rbp_kd_lib_list[grep('HepG2', rbp_kd_lib_list)]]
k562_minus_library <- kd_library_all[names(kd_library_all) != rbp_kd_lib_list[grep('K562', rbp_kd_lib_list)]]

# add new study
hepg2_minus_library[[input_study]] <- study_mats
k562_minus_library[[input_study]] <- study_mats

## run analysis
hepg2_study <- rbp_kd_lib_list[grep('HepG2', rbp_kd_lib_list)]
hepg2_test_result_df <- perform_SPARKS_analysis_with_overlap_filter(kd_library_all[[hepg2_study]],
                                                                    hepg2_minus_library,
                                                                    study = hepg2_study,
                                                                    num_cores = 4,
                                                                    overlap_ratio_threshold = num_ratio)
hepg2_output_table <- sprintf("./swap/benchmark_result.%s.%s.%s.HepG2_swapped.txt", input_study, spl_type, num_events)

data.table::fwrite(hepg2_test_result_df,
                   file = hepg2_output_table)

# k562 run
k562_study <- rbp_kd_lib_list[grep('K562', rbp_kd_lib_list)]
k562_test_result_df <- perform_SPARKS_analysis_with_overlap_filter(kd_library_all[[k562_study]],
                                                                    k562_minus_library,
                                                                    study = k562_study,
                                                                    num_cores = 4,
                                                                    overlap_ratio_threshold = num_ratio)
k562_output_table <- sprintf("./swap/benchmark_result.%s.%s.%s.K562_swapped.txt", input_study, spl_type, num_events)

data.table::fwrite(k562_test_result_df,
                   file = k562_output_table)

# run normal one
# k562 run
original_test_result_df <- perform_SPARKS_analysis_with_overlap_filter(study_mats,
                                                                   kd_library_all,
                                                                   study = input_study,
                                                                   num_cores = 4,
                                                                   overlap_ratio_threshold = num_ratio)
original_output_table <- sprintf("./swap/benchmark_result.%s.%s.%s.original.txt", input_study, spl_type, num_events)

data.table::fwrite(original_test_result_df,
                   file = original_output_table)
