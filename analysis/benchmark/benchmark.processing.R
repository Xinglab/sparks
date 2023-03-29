##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(fgsea)
library(metap)
library(dplyr)


##### ANALYSIS CODE #####
# import library for analysis
print("Loading Database")

spl_types <- c("SE", "A3SS", "A5SS")
# spl_type  <- "A5SS"
spl_type <- "SE"

kd_library_all <- readRDS(sprintf('/home/yangt3/xinglab/SPARKS/library/Clean_minfiltered.%s.library.rds', spl_type))

output_dir <- "./"
count_threshold <- 20

##### ANALYSIS #####
## stratify the experiments based on the method
# perturbation_methods <- c("shRNA", "CRISPR")
perturbation_method <- "shRNA"

# all_study_list <- names(kd_library_all[["SE"]])
all_study_list <- names(kd_library_all)


print(sprintf("Processing %s", perturbation_method))

# subset the library by name
input_study_list <- all_study_list[grep(perturbation_method, all_study_list)]

## identify replicated RBPs in the two cell lines
# generate study -> RBP table
rbp_info_df <- data.frame(study = input_study_list,
                          rbp = unlist(lapply(input_study_list, function(x) strsplit(x, "_")[[1]][2])))
rbp_num_exp <- table(rbp_info_df$rbp)
rbp_kd2_list <- names(rbp_num_exp)[rbp_num_exp == 2]

encode_rep_experiments <- rbp_info_df[rbp_info_df$rbp %in% rbp_kd2_list, ]$study


### Perform Benchmark
# input the argument
args = commandArgs(trailingOnly=TRUE)
study_of_interest <- args[1]
num_events <- as.integer(args[2])

num_ratio <- num_events / 100

# specify the data for the splicing type
study_cell_line <- strsplit(study_of_interest, "_")[[1]][1]

other_cell_line_experiments <- encode_rep_experiments[!(startsWith(encode_rep_experiments, study_cell_line))]

### Perform Benchmark

print(study_of_interest)

# define the test dataset
study_mats <- kd_library_all[[study_of_interest]]

# define the output file
# - we are only saving the results since saving object would spend too much storage
output_table <- sprintf("../result/benchmark_result.new_GSEA_rev.%s.%s.%s.filtered_fgsea.%s.df.txt", perturbation_method, spl_type, study_of_interest, num_events)

# run analysis if the file is already not there
if (!(file.exists(output_table))){

  # run analysis
  # - SPARKS analysis code actually include all three similarity metrics
  # - Note that the rank would be wrong, so would need to re compute downstream
  test_result_df <- perform_SPARKS_analysis_with_overlap_filter(study_mats,
                                                                kd_library_all,
                                                                study = study_of_interest,
                                                                num_cores = 4,
                                                                overlap_ratio_threshold = num_ratio,
                                                                library_list = encode_rep_experiments)
  print(dim(test_result_df))
  data.table::fwrite(test_result_df,
                     file = output_table)

}
