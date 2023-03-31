##### FUNCTIONS #####
#' convert MATS to SPARKS format
#'
#' @param input_mats filtered ${spl_type}.MATS.JC.txt
#'
#' @return converted MATS
#' @export
convert_MATS_format <- function(input_mats){
  new_mats <- do.call(rbind, apply(input_mats, 1, function(event_entry){
    # combine event ID
    event_id <- paste(c(trimws(event_entry[2:11]), spl_type), collapse = ":")

    # remove unnecessary coordinates
    event_id_clean <- SPARKS::rewrite_event_coordinates(event_id)

    # parse comma separated count info
    inc_list_group1 <- as.numeric(strsplit(as.character(event_entry['IJC_SAMPLE_1']), split = ',')[[1]])
    skip_list_group1 <- as.numeric(strsplit(as.character(event_entry['SJC_SAMPLE_1']), split = ',')[[1]])

    inc_list_group2 <- as.numeric(strsplit(as.character(event_entry['IJC_SAMPLE_2']), split = ',')[[1]])
    skip_list_group2 <- as.numeric(strsplit(as.character(event_entry['SJC_SAMPLE_2']), split = ',')[[1]])

    # pool the counts in each replicates
    inc_count_group1 <- sum(inc_list_group1)
    skip_count_group1 <- sum(skip_list_group1)

    inc_count_group2 <- sum(inc_list_group2)
    skip_count_group2 <- sum(skip_list_group2)

    # query length for psi calculation
    inc_len <- as.numeric(event_entry['IncFormLen'])
    skip_len <- as.numeric(event_entry['SkipFormLen'])

    # calcaulte new PSI
    psi_group1 = (inc_count_group1 / inc_len) / ((inc_count_group1 / inc_len) + (skip_count_group1 / skip_len))
    psi_group2 = (inc_count_group2 / inc_len) / ((inc_count_group2 / inc_len) + (skip_count_group2 / skip_len))

    # calculate delta PSI
    delta_psi_pulled = psi_group1 - psi_group2

    # calculate total, avg, min count
    inc_count <- c(inc_list_group1, inc_list_group2)
    skip_count <- c(skip_list_group1, skip_list_group2)

    total_count <- inc_count + skip_count
    avg_count <- mean(total_count)
    min_count <- min(total_count)

    # merge raw counts for downstream analysis
    count_values <- paste0(total_count, collapse = ",")

    # merge raw PSI values
    psi_values <- paste0(c(event_entry["IncLevel1"], event_entry["IncLevel2"]), collapse = ",")

    # extract relevant information
    # beta <- event_entry['IncLevelDifference']
    beta <- delta_psi_pulled
    pval <- event_entry['PValue']
    fdr <- event_entry['FDR']

    # combine the event info as a data
    clean_mats_entry <- data.frame(event_id_clean,
                                   as.numeric(beta),
                                   as.numeric(pval),
                                   as.numeric(fdr),
                                   as.numeric(avg_count),
                                   psi_values,
                                   count_values,
                                   as.numeric(inc_len),
                                   as.numeric(skip_len),
                                   as.numeric(psi_group1),
                                   as.numeric(psi_group2),
                                   as.numeric(delta_psi_pulled),
                                   as.numeric(min_count))

    colnames(clean_mats_entry) <- c('event',
                                    'beta',
                                    'pval',
                                    'fdr',
                                    'avg_count',
                                    'psi_values',
                                    'count_values',
                                    'inc_len',
                                    'skip_len',
                                    'pulled_psi_1',
                                    'pulled_psi_2',
                                    'pulled_delta_psi',
                                    'min_count')

    return(clean_mats_entry)
  }))
  return(new_mats)
}

#' Import raw rMATS output
#'
#' Process raw rMATS output for SPARKS
#'
#' @param input_dir directory with rMATS Output
#' @param spl_type splicing type (default = SE)
#'
#' @return mats_df
#' @export
import_raw_rMATS_output <- function(input_dir, spl_type = "SE"){
  # read in the mats
  input_file <- sprintf("%s/%s.MATS.JC.txt", input_dir, spl_type)
  input_mats <- data.table::fread(input_file)

  # read in fromGTF files to determine known JC files
  input_event_file_all <- sprintf("%s/fromGTF.%s.txt", input_dir, spl_type)
  input_event_file_novel_jxn <- sprintf("%s/fromGTF.novelJunction.%s.txt", input_dir, spl_type)
  input_event_file_novel_ss <- sprintf("%s/fromGTF.novelSpliceSite.%s.txt", input_dir, spl_type)

  input_event_all <- data.table::fread(input_event_file_all)$ID
  input_event_novel_jxn <- data.table::fread(input_event_file_novel_jxn)$ID
  input_event_novel_ss <- data.table::fread(input_event_file_novel_ss)$ID

  # remove all IDs found on novel jxn or novel ss from all events
  known_ids <- input_event_all[!(input_event_all %in% union(input_event_novel_jxn,
                                                            input_event_novel_ss))]

  # filter mats by ID
  filtered_mats <- input_mats[input_mats$ID %in% known_ids, ]

  # convert the format to SPARKS-usable one
  new_mats <- convert_MATS_format(filtered_mats)

  # sort and filter
  min_filtered_mats <- subset(new_mats, min_count >= 20)

  sorted_mats <- min_filtered_mats %>% arrange(-beta)

  return(sorted_mats)
}

