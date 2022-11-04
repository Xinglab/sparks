library(START)
# TODO - generalize it for other splicing types
spl_type <- "SE"

# find all the input files
# local debug
# input_dir <- "/Users/harryyang/research/Xing_Lab/START/analysis/CLIP"
# clip_exp_list <- read.csv('/Users/harryyang/research/Xing_Lab/START/analysis/clip_exp_list.txt', header = F)$V1
## input_dir <- "/home/yangt3/xinglab/ENCODE/START/run/library/" # TODO
## clip_exp_list <- read.csv('/home/yangt3/xinglab/ENCODE/bed/clip_exp_list.txt', header = F)$V1 # TODO
input_start_files <- Sys.glob(sprintf("%s/*.START.rds", input_dir))

# make lists for processings
event_list <- c()

region_list <- c("Upstream_5ss_exon", "Upstream_5ss_intron",
                 "Cassette_3ss_intron", "Dnstream_3ss_exon",
                 "Dnstream_3ss_intron", "Cassette_5ss_intron",
                 "Cassette_3ss_exon", "Cassette_5ss_exon")
region_clip_list <- list()
dummy <- lapply(region_list, function(region){
  region_clip_list[[region]] <<- list()
  return(NULL)
})

dummy <- lapply(input_start_files, function(input_start_file){
  print(input_start_file)
  input_start <- readRDS(input_start_file)

  # extract new events not covered by previous ones
  # - as the peaks are coming from coordinates, it will be the same regardless of the experiment
  events <- rownames(input_start@clip_result_list[[spl_type]]$Upstream_5ss_exon)
  known_events <- rownames(subset(input_start@exon_annotation[[spl_type]], annotation == "Known_JC"))

  clip_events <- events[events %in% known_events]
  new_events <- clip_events[!(clip_events %in% event_list)]

  # add the new events to the list
  event_list <<- c(event_list, new_events)

  dummy <- lapply(region_list, function(region){
    clip_result <- input_start@clip_result_list[[spl_type]][[region]]

    # extract new events
    clip_result_select <- clip_result[new_events, ]

    # add null results - if no CLIP peak is found, these are omitted,
    # thus it is necessary to add them for downstream merge,
    # although there is no peak for this particular KD experiment
    null_clip_list <- clip_exp_list[!(clip_exp_list %in% colnames(clip_result_select))]

    clip_result_select[, null_clip_list] <- "no_peak"

    # add the new ones to the list
    region_clip_list[[region]][[length(region_clip_list[[region]]) + 1]] <<- clip_result_select
    return(NULL)
  })

  return(NULL)
})

compiled_region_list <- list()
dummy <- lapply(region_list, function(region){

  region_clip_df <- do.call(rbind, region_clip_list[[region]])

  compiled_region_list[[region]] <<- region_clip_df
  return(NULL)
})

saveRDS(compiled_region_list, sprintf("./CLIP_result_library.%s.rds", spl_type))
