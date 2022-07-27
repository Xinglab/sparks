##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)

# import library for analysis
kd_library_all <- readRDS('/home/yangt3/xinglab/ENCODE/SPARKS/ENCODE_library.KD_and_KO.count_20.rds')

output_dir <- "./"
count_threshold <- 20

spl_types <- c("SE", "A3SS", "A5SS")

##### ANALYSIS #####
## stratify the experiments based on the method
perturbation_methods <- c("shRNA", "CRISPR")

all_study_list <- names(kd_library_all[["SE"]])

### Perform Benchmark
dummy <- lapply(spl_types, function(spl_type) {
  print(spl_type)

  # specify the data for the splicing type
  kd_library_spltype <- kd_library_all[[spl_type]]

  benchmark_result_df <- do.call(rbind, lapply(all_study_list, function(study_of_interest){
    # define the test dataset
    study_mats <- kd_library_spltype[[study_of_interest]]
    print(study_of_interest)

    # run analysis
    # - SPARKS analysis code actually include all three similarity metrics
    # - Note that the rank would be wrong, so would need to re compute downstream
    test_result_df <- SPARKS::perform_SPARKS_analysis(study_mats, kd_library_spltype, study = study_of_interest, num_cores = 24)

    return(test_result_df)
  }))

  data.table::fwrite(benchmark_result_df,
                     file = sprintf("./benchmark/ENCODE.pairwise_all.%s.df.txt", spl_type))
  return(count_threshold)
})

