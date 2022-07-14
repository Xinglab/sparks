##### INITIALIZATION #####
library(START)
library(ggplot2)
library(pbmcapply)

output_dir <- "./"
count_threshold = 20

import_START_MATS_for_analysis <- function(input_start, event_type, count_threshold){
  # process the data for analysis
  known_events <- rownames(subset(input_start@exon_annotation[[event_type]], annotation == "Known_JC"))
  study_mats_raw <- input_start@MATS_list[[event_type]]

  study_mats_known <- study_mats_raw[study_mats_raw$event %in% known_events, ]
  study_mats_temp <- subset(study_mats_known, avg_count >= count_threshold)
  return(study_mats_temp)
}

### read in the data

## input data list
# input_start_dir <- "/Users/harryyang/research/Xing_Lab/START/test"
input_start_dir <- "/home/yangt3/xinglab/ENCODE/START/run/library/"

input_start_list <- Sys.glob(sprintf('%s/*.START.rds', input_start_dir))

# get study list from file names
input_study_list <- unlist(lapply(input_start_list, function(input_file){
  # define name
  input_file_name <- strsplit(input_file, split = "/")[[1]][length(strsplit(input_file, split = "/")[[1]])]
  input_study <- strsplit(input_file_name, "\\.")[[1]][1]

  return(input_study)
}))
# rename for easier lookup
names(input_start_list) <- input_study_list


# apply filter to all KD mats
# discordant_mats_list <- list()
concordant_mats_list <- list()

# make slots
concordant_mats_list[['SE']] <- list()
concordant_mats_list[['A3SS']] <- list()
concordant_mats_list[['A5SS']] <- list()
concordant_mats_list[['RI']] <- list()

spl_type_list <- c("SE", "A3SS", "A5SS", "RI")
## import START file for processing
# these are sorted by beta
dummy <- lapply(input_study_list, function(study){
  print(study)
  input_start_file <- input_start_list[study]

  # read both files
  input_start <- readRDS(input_start_file)



  for (spl_type in spl_type_list){
    print(spl_type)
    # input_mats <- input_start@MATS_list[[spl_type]]


    # input_mats$min_count <- unlist(lapply(input_mats$count_values, function(x) min(do.call(as.numeric, strsplit(x, ",")))))
    # input_mats_filtered <- subset(input_mats, min_count >= count_threshold)

    input_mats_filtered <- import_START_MATS_for_analysis(input_start, event_type = spl_type, count_threshold = count_threshold)

    input_mats_sorted <- input_mats_filtered[rev(order(input_mats_filtered$beta)), ]

    # ggplot(input_mats_dump, aes(x = log2(avg_count + 1))) + geom_density() + geom_vline(xintercept = log2(20 + 1)) + facet_wrap(~ifelse(abs(pulled_delta_psi) >= 0.05, "Sig dif", "nonsig"), scale = "free_y")

    # segregate discordant vs. concordant list
    # discordant_mats <- subset(input_mats_filtered, abs(beta - pulled_delta_psi) > 0.1)
    concordant_mats <- subset(input_mats_sorted, abs(beta - pulled_delta_psi) < 0.1)

    # change beta
    concordant_mats$mean_beta <- concordant_mats$beta
    # discordant_mats$mean_beta <- discordant_mats$beta
    concordant_mats$beta <- concordant_mats$pulled_delta_psi
    # discordant_mats$beta <- discordant_mats$pulled_delta_psi

    # sort by beta
    # concordant_mats_sorted <- concordant_mats[rev(order(concordant_mats$beta)), ]

    concordant_mats_list[[spl_type]][[study]] <<- concordant_mats
    # discordant_mats_list[[study]] <<- discordant_mats

  }
  return()
})

# save the signature
output_file <- sprintf("%s/ENCODE.shRNA_KD.signatures.MATS_df_list.count_%s.rds", output_dir, count_threshold)

saveRDS(concordant_mats_list, output_file)


