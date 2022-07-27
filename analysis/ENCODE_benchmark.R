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

lapply(perturbation_methods, function(perturbation_method){
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

  print(encode_rep_experiments)

  ### Perform Benchmark
  dummy <- lapply(spl_types, function(spl_type) {
    print(spl_type)

    # specify the data for the splicing type
    kd_library_spltype <- kd_library_all[[spl_type]][encode_rep_experiments]

    benchmark_result_df <- do.call(rbind, lapply(encode_rep_experiments, function(study_of_interest){
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
                       file = sprintf("./benchmark/benchmark_result.%s.%s.df.txt", perturbation_method, spl_type))
    return(count_threshold)
  })

})
