usethis::use_package("R.utils")
usethis::use_package("data.table")
usethis::use_package("dplyr")
usethis::use_package("reshape2")
usethis::use_package("gtools")
usethis::use_package("maftools")
usethis::use_package("stringr")
usethis::use_package("ggplot2")
usethis::use_package("RColorBrewer")
usethis::use_package("pbmcapply")

#' @export
SPARKS <- setClass("SPARKS",
                  slots = c(
                    study = "character",
                    num_sample = "numeric",
                    num_mutation = "numeric",
                    sample_chart = "data.frame",
                    # below are list - each list would have slots for each splicing event
                    psi_df = "list",
                    linear_model_pval = "list",
                    linear_model_beta = "list",
                    motif_enrichment_df = "list",
                    # below are dataframe as these are common among the splicing events
                    hotspot_maf = "MAF",
                    hotspot_mutation_matrix = "data.frame",
                    exp_df = "data.frame",
                    # TODO - add DESeq part,
                    diff_exp_result = "data.frame",
                    deseq_result = "data.frame",
                    # below are parts of data from analysis
                    summary_asso_df = "data.frame",
                    motif_analysis_result = "list",
                    clip_result_list = "list",
                    exon_annotation = "list",
                    MATS_list = "list"
                  ))


##### SUBFUNCTIONS #####


#' Clean up pval data for downstream analysis
#'
#' This function nicely gather and sort hotspot and AS events by coordinates
#' then convert the information into neglog10 for easier handling
#'
#'
#' @param input_pval_df pval dataframe form analysis script
#' @param spl_event_type either SE, A3SS, A5SS
#' @return pval_df where rows are AS events, columns are mutation/sample, elements are -log10(pval)
#' @export
cleanup_hotspot_pval_data <- function(input_pval_df, spl_event_type){
  ### This function cleans up data
  ## nicely gather and sort hotspot and AS events by coordinates
  ## then convert the information into neglog10 for easier handling
  # spl_event_type can be either SE, A3SS, A5SS
  # spl event type will be added as additional column


  ## de-parse column labels - first 10 rows
  # columns are AS event
  # first 6 columns are excluded as these are for hotspot info and empty for AS event
  pval_df_cols <- as.data.frame(t(input_pval_df[seq(10), -seq(6)]))

  # add splicing type for AS labels
  pval_df_cols$spl <- spl_event_type

  # concat using ":"
  pval_df_AS_labels <- as.vector(apply(pval_df_cols, 1, function(x) paste0(x, collapse = ":")))

  ## generate order for heatmap sorting later
  coords <- do.call(cbind, lapply(seq(5, 10), function(x) as.numeric(pval_df_cols[, x])))
  mean_coord <- rowMeans(coords)

  # bind chromosome and mean coordinate
  chr_coord <- apply(cbind(pval_df_cols[, 3], mean_coord), 1,
                     function(x) paste0(x, collapse = ":"))
  names(chr_coord) <- pval_df_AS_labels

  # mixsort to find the sorted names
  chr_coord_sorted <- gtools::mixedsort(chr_coord)


  ## de-parse row labels - first 6 columns
  # rows are hotspot
  pval_df_rows <- input_pval_df[-seq(10), seq(6)]

  # reorganize column orders
  colnames(pval_df_rows) <- as.vector(t(pval_df_rows[1,]))
  pval_df_rows <- pval_df_rows[-1, ]

  pval_df_rows <- pval_df_rows[, c("geneID", "geneSymbol", "chromosome",
                                   "hotspot_start_pos", "hotspot_end_pos",
                                   "count")]

  # concat using ":"
  pval_df_HS_labels <- unlist(as.vector(apply(pval_df_rows, 1, function(x) paste0(x, collapse = ":"))))

  ## find order
  hs_chr_coord <- apply(pval_df_rows[, c('chromosome', 'hotspot_start_pos')],
                        1, function(x) paste0(x, collapse = ":"))

  names(hs_chr_coord) <- pval_df_HS_labels

  hs_chr_coord_sorted <- mixedsort(hs_chr_coord)
  # extract data part
  pval_df <- data.frame(input_pval_df[-seq(11), -seq(6)], stringsAsFactors = FALSE)
  colnames(pval_df) <- pval_df_AS_labels
  rownames(pval_df) <- pval_df_HS_labels

  pval_df <- pval_df[names(hs_chr_coord_sorted), names(chr_coord_sorted)]

  # convert the p values into numbers
  pval_df[] <- lapply(pval_df, as.character)
  pval_df[] <- lapply(pval_df, as.numeric)

  return(pval_df)
}




#' Remove duplicates arising from sample processing
#'
#' This function removes columns (samples) in PSI df where there may be duplicates
#' generated from some error in processing pipeline
#'
#' @param psi_df PSI value dataframe
#' @param event_type either SE, A3SS, A5SS
#' @return clean_psi_df - PSI df without duplicated samples
#' @export
cleanup_psi_df <- function(psi_df, event_type){
  psi_rownames <- unlist(apply(psi_df[, seq(10)], 1, function(x) gsub(paste0(c(x, event_type), collapse = ":"), pattern = " ", replacement = "")))

  # clean up row names
  psi_df_filtered <- psi_df[,-seq(10)]
  rownames(psi_df_filtered) <- psi_rownames

  # drop dup columns
  psi_df_selected <- psi_df_filtered[, !unlist(lapply(colnames(psi_df_filtered), function(x) endsWith(x, "A.")))]

  return(psi_df_selected)
}




#' Remove duplicates arising from sample processing
#'
#' This function removes columns (samples) in EXP df where there may be duplicates
#' generated from some error in processing pipeline
#'
#' @param exp_df EXP value dataframe
#' @return clean_exp_df - PSI df without duplicated samples
#' @export
cleanup_exp_df <- function(psi_df){
  event_type <- "Expression"
  psi_rownames <- unlist(apply(psi_df[, seq(2)], 1, function(x) gsub(paste0(c(x, event_type), collapse = ":"), pattern = " ", replacement = "")))

  # clean up row names
  psi_df_filtered <- psi_df[,-seq(2)]
  rownames(psi_df_filtered) <- psi_rownames

  # drop dup columns
  psi_df_selected <- psi_df_filtered[, !unlist(lapply(colnames(psi_df_filtered), function(x) endsWith(x, "A.")))]

  return(psi_df_selected)
}



#' Wrapper for cleaning up hotspot pval dataframe
#'
#' This function reads in the file and cleans up accordingly
#'
#' @param input_hotspot_lm_file hotspot lm result file from the processing pipeline
#' @param spl_event_type either SE, A3SS, A5SS
#' @return clean_pval_df - pval df where col/rownames are processed
#' @export
import_linear_model_raw_data <- function(input_hotspot_lm_file,
                                         spl_event_type){
  input_pval_df <- read.csv(input_hotspot_lm_file, header = FALSE)
  # clean up annotation columns and rows
  clean_pval_df <- cleanup_hotspot_pval_data(input_pval_df, spl_event_type = spl_event_type)

  # associations will be jointly filtered with beta and pval
  # - this code can be used for both beta and pval import
  return(clean_pval_df)
}




#' Wrapper for cleaning up hotspot mutation type
#'
#' This function parses hotspot mutations for downstream visualization
#'
#' @param hs_mut_all hotspot mutation
#' @return sorted hotspot mutation matrix
#' @export
cleanup_hotspot_mutation_matrix <- function(hs_mut_all){
  # import hotspot mutation matrix for downstream analysis

  # extract relavent columns for labeling
  hs_rows <- hs_mut_all[, seq(6)]

  # reorganize column orders
  hs_rows <- hs_rows[, c("geneID", "geneSymbol", "chromosome",
                         "hotspot_start_pos", "hotspot_end_pos", "count")]

  # concat using ":"
  # Whitespace needs to be trimmed for data cleanup
  hs_rows_label <- unlist(as.vector(apply(hs_rows, 1, function(x) paste0(unlist(lapply(x, function(y) trimws(y))), collapse = ":"))))

  # extract hotspot-data matrix
  hs_mut_df <- hs_mut_all[, -seq(6)]
  rownames(hs_mut_df) <- hs_rows_label

  ## find order
  # Whitespace needs to be trimmed for correct sorting
  hs_chr_coord <- apply(hs_rows[, c('chromosome', 'hotspot_start_pos')],
                        1, function(x) paste0(unlist(lapply(x, function(y) trimws(y))), collapse = ":"))

  names(hs_chr_coord) <- hs_rows_label

  # sort them by genomic positions
  hs_chr_coord_sorted <- mixedsort(hs_chr_coord)

  # sort the dataframes by genomic coordinates
  hs_mut_df <- hs_mut_df[names(hs_chr_coord_sorted), ]

  # fill NA with 0
  hs_mut_df[is.na(hs_mut_df)] <- 0

  return(hs_mut_df)

}



#' Internal update function to update the sample chart
#'
#' This function updates the sample chart to add extra annotations arising from analysis
#'
#' @param object the object
#' @param data_name column name to annotate the data
#' @param samples samples that have certain features - 0/1 will be added accordingly
#' @return updated object with additional annotation
#' @export
update_sample_chart <- function(object, data_name, samples){
  if (dim(object@sample_chart)[1] != 0){
    # cross-reference sample names
    sample_indicator <- as.data.frame(ifelse(rownames(object@sample_chart) %in% samples, 1, 0))
    colnames(sample_indicator) <- c(data_name)

    # bind the column with sample chart
    new_sample_chart <- cbind(object@sample_chart, sample_indicator)
    object@sample_chart <- new_sample_chart
    # TODO - this method does not account for cases where there is sample that
    # are not in the sample chart are present. This needs to be fixed.

  } else { # if there is no sample chart, it needs to be contructed
    sample_indicator <- as.data.frame(rep(1, length(samples)))
    rownames(sample_indicator) <- samples
    colnames(sample_indicator) <- c(data_name)

    object@sample_chart <- sample_indicator
  }
  return(object)
}



#' Parse COSMIC data
#'
#' Parses COSMIC genes to classify either Tumor suppressing or oncogene
#'
#' @param cosmoc_df COSMIC df read in from COSMIC data from website
#' @return list of genes grouped in TSG/OG/None separately
#' @export
classify_cosmic_genes <- function(cosmic_df){

  # prepare the oncogene/tumor suppressor gene
  cosmic_df$TSG <- ifelse(grepl("TSG", cosmic_df$Role.in.Cancer), 1, 0)
  cosmic_df$OG <- ifelse(grepl("oncogene", cosmic_df$Role.in.Cancer), 1, 0)

  ## resolve aliases for genes
  # extract synonymous TSG and OG
  cosmic_list_tsg <- unlist(lapply(cosmic_df[cosmic_df$TSG == 1, ]$Synonyms, function(x) strsplit(split = ",", as.character(x))))
  cosmic_list_og <- unlist(lapply(cosmic_df[cosmic_df$OG == 1, ]$Synonyms, function(x) strsplit(split = ",", as.character(x))))

  # add the original genes
  cosmic_list_tsg <- c(cosmic_df[cosmic_df$TSG == 1, ]$Gene.Symbol,
                       cosmic_list_tsg)
  cosmic_list_og <- c(cosmic_df[cosmic_df$OG == 1, ]$Gene.Symbol,
                      cosmic_list_og)

  # clean up by removing duplicates
  cosmic_list_tsg <- unique(cosmic_list_tsg)
  cosmic_list_og <- unique(cosmic_list_og)

  # find genes that are in neither tsg or og - probably fusion
  cosmic_list_none <- as.character(cosmic_df$Gene.Symbol[!(cosmic_df$Gene.Symbol
                                                           %in% cosmic_list_tsg |
                                                             cosmic_df$Gene.Symbol
                                                           %in% cosmic_list_og)])

  cosmic_list <- list("TSG" = cosmic_list_tsg,
                      "OG" = cosmic_list_og,
                      "None" = cosmic_list_none)

  return(cosmic_list)
}




##### MAIN IMPORT FUNCTIONS #####


#' Main import function to import PSI df into object
#'
#' This function reads in the file and cleans up accordingly and add it to the object
#'
#' @param object the object
#' @param input_psi_file psi file from the processing pipeline
#' @param spl_event_type either SE, A3SS, A5SS
#' @return object with the PSI dataframe added
#' @export
import_PSI_df <- function(object, input_psi_file, spl_event_type){
  print(sprintf("Importing PSI data for %s", spl_event_type))
  input_df <- read.csv(input_psi_file, check.names = F, stringsAsFactors = F)

  # clean up the data for usable input format
  clean_df <- cleanup_psi_df(input_df, spl_event_type)

  # add the clean PSI df to object
  object@psi_df[[spl_event_type]] <- clean_df

  # update the sample chart to keep track of samples
  samples <- colnames(clean_df)
  object <- update_sample_chart(object, sprintf("PSI_%s", spl_event_type), samples)
  return(object)
}



#' Main import function to import pval df into object
#'
#' This function reads in the file and filteres the data and add it to the object
#'
#' @param object the object
#' @param input_hotspot_lm_pval_file pval file from the processing pipeline
#' @param input_hotspot_lm_beta_file beta file from the processing pipeline
#' @param spl_event_type either SE, A3SS, A5SS
#' @return object with the linear model result dataframe added
#' @export
import_linear_model_data <- function(object,
                                     input_hotspot_lm_pval_file,
                                     input_hotspot_lm_beta_file,
                                     spl_event_type,
                                     neglog_pval_cutoff = 5,
                                     beta_cutoff = 0.1){
  # import files
  print(sprintf("Importing Linear Model results for %s", spl_event_type))
  pval_df <- import_linear_model_raw_data(input_hotspot_lm_pval_file, spl_event_type)
  beta_df <- import_linear_model_raw_data(input_hotspot_lm_beta_file, spl_event_type)

  # validate the data
  if (sum(colnames(pval_df) == colnames(beta_df))/length(colnames(pval_df)) != 1){
    stop("Splicing events in beta and pval data do not match - check your input")
  }
  if (sum(rownames(pval_df) == rownames(beta_df))/length(rownames(pval_df)) != 1){
    stop("Mutation events in beta and pval data do not match - check your input")
  }

  ## filter the events
  # create mask
  pval_mask <- pval_df
  pval_mask[pval_mask < neglog_pval_cutoff] <- 0
  pval_mask[pval_mask > neglog_pval_cutoff] <- 1

  beta_mask <- beta_df
  beta_mask[abs(beta_mask) < beta_cutoff] <- 0
  beta_mask[abs(beta_mask) > beta_cutoff] <- 1

  # combine mask to select events that pass both filter
  combined_mask <- pval_mask + beta_mask

  # flatten combined mask to use for multiplication
  combined_mask[combined_mask < 2] <- 0
  combined_mask[combined_mask >= 2] <- 1

  # filter pval and beta by the mask
  pval_df_masked <- combined_mask * pval_df
  beta_df_masked <- combined_mask * beta_df

  # rmeove unnecessary rows in both pval and beta df with 0 informative row/cols
  pval_df_filtered <- pval_df_masked[rowSums(pval_df_masked) != 0, colSums(pval_df_masked) != 0]
  beta_df_filtered <- beta_df_masked[rowSums(beta_df_masked) != 0, colSums(beta_df_masked) != 0]

  # save those information to the object
  object@linear_model_pval[[spl_event_type]] <- pval_df_filtered
  object@linear_model_beta[[spl_event_type]] <- beta_df_filtered

  return(object)
}



#' Main import function to import MAF into object
#'
#' This function reads in the MAF file and add it to the object
#'
#' @param object the object
#' @param input_maf_file maf file from analysis pipeline - hotspot should be added in this MAF
#' @return object with the MAF result dataframe added
#' @export
import_hotspot_maf <- function(object,
                               input_maf_file){


  # MAF object from maftools is used for downstream plotting purposes like lollipop plot
  input_maf <- read.maf(input_maf_file)

  object@hotspot_maf <- input_maf
  return(object)
}



#' Main import function to import expression df into object
#'
#' This function reads in the file and filteres the data and add it to the object
#'
#' @param object the object
#' @param input_exp_file expression matrix file from the processing pipeline
#' @return object with the expression dataframe added
#' @export
import_expression_matrix <- function(object,
                                     input_exp_file){
  print("Importing Exp data")
  input_df <- read.csv(input_exp_file, check.names = F, stringsAsFactors = F, sep='\t')

  # clean up the data for usable input format
  clean_df <- cleanup_exp_df(input_df)
  # add the clean PSI df to object
  object@exp_df <- clean_df

  # add annotation to sample chart
  samples <- colnames(clean_df)
  object <- update_sample_chart(object, "EXP", samples)
  return(object)
}



#' Main import function to import hotspot matrix into object
#'
#' This function reads in the file and add it to the object
#' Hotspot matrix is where columns are sample, rows are hotspot mutations,
#' 1 if the sample has mutation in the hotspot, 0 if not
#' @param object the object
#' @param input_hotspot_matrix_file hotspot matrix file from the processing pipeline
#' @return object with the expression dataframe added
#' @export
import_hotspot_mutation_matrix <- function(object,
                                           input_hotspot_matrix_file){
  print("Importing Hotspot Mutation Matrix")
  # read in the files
  input_hotspot_df <- read.csv(input_hotspot_mutation_file, check.names = F)
  clean_hotspot_df <- cleanup_hotspot_mutation_matrix(input_hotspot_df)

  object@hotspot_mutation_matrix <- clean_hotspot_df

  # update the sample chart
  samples <- colnames(clean_hotspot_df)
  object <- update_sample_chart(object, "Hotspot_Mutation", samples)

  return(object)
}




#' Main import function to import motif count into object
#'
#' This function reads in the file and add it to the object

#' @param object the object
#' @param motif_count_dir directory with motif count data generated putative RBP binding motifs
#' @param spl_event_type event_type
#' @return object with the list of dataframe added
#' @export
import_motif_count_data <- function(object, motif_count_dir, spl_event_type){
  input_peak_file_list <- Sys.glob(sprintf("%s/score.*.%s.txt", motif_count_dir, spl_event_type))

  region_count_list = list()
  dummy <- lapply(input_peak_file_list, function(input_peak_file){
    peak_df <- as.data.frame(fread(input_peak_file))

    # extract motif name
    peak_file_split <- strsplit(input_peak_file, '/')[[1]]
    motif_name <- strsplit(peak_file_split[length(peak_file_split)], '\\.')[[1]][2]

    # progress report
    print(sprintf("Importing %s on %s", motif_name, spl_event_type))

    # prepare the data
    rownames(peak_df) <- peak_df$ID
    peak_df$ID <- NULL

    # split the dataframe into each region
    for (region in colnames(peak_df)[1:length(colnames(peak_df))-1]){

      region_count_df = as.data.frame(peak_df[, region])
      rownames(region_count_df) <- rownames(peak_df)
      colnames(region_count_df) <- motif_name

      # add new entries to append to if the region_count_list is empty
      if (is.null(region_count_list[[region]])){
        region_count_list[[region]] = list()
      }

      # add the entry
      region_count_list[[region]][[motif_name]] <<- region_count_df
    }
  })

  # merge the dataframe
  count_matrix_list <- list()
  dummy <- lapply(names(region_count_list), function(region){
    region_peak_count_matrix <- as.data.frame(do.call(cbind, region_count_list[[region]]))
    count_matrix_list[[region]] <<- region_peak_count_matrix
  })

  # assign the count matrix to the object
  object@motif_enrichment_df[[spl_event_type]] <- count_matrix_list

  return(object)
}




#' Main import function to import ENCODE CLIP peak count into object
#'
#' This function reads in the file and add it to the object

#' @param object the object
#' @param clip_result_file CLIP peak overlap file generated from analysis pipeline
#' @param event_type event_type
#' @return object with the list of dataframe added
#' @export
import_ENCODE_CLIP_intersect_data <- function(object, clip_result_file, event_type){
  clip_result <- fread(clip_result_file, header = TRUE, sep='\t')
  region_list <- unique(clip_result$region)
  asso_df_region_list <- lapply(region_list, function(region_of_interest){

    asso_df_list <- pbmclapply(as.data.frame(t(unique(clip_result[,c("study", "RBP")]))), function(entry){
      study_of_interest <- entry[1]
      rbp_of_interest <- entry[2]
      clip_subset = subset(clip_result, study == study_of_interest &
                             RBP == rbp_of_interest &
                             region == region_of_interest)

      asso_df <- as.data.frame(rownames(object@psi_df[[event_type]])) # retreive all events - serves as background
      colnames(asso_df) <- c("AS_events")
      asso_df[[sprintf("%s_%s", study_of_interest, rbp_of_interest)]] <- ifelse(asso_df$AS_events %in% clip_subset$exon, "peak", "no_peak")

      rownames(asso_df) <- asso_df$AS_events
      asso_df$AS_events <- NULL
      return(asso_df)
    }, mc.cores = 4)

    clip_peak_df <- do.call(cbind, asso_df_list)
    return(clip_peak_df)
  })
  names(asso_df_region_list) <- region_list

  # assign the result to object
  object@clip_result_list[[event_type]] <- asso_df_region_list
  return(object)
}





#' Main import function to import exon annotation into object
#'
#' This function reads in the file and add it to the object

#' @param object the object
#' @param annotation_file Exon annotation file (KnownJC/NovelJC/NovelSS) annotated for each AS event
#' @param event_type event_type
#' @return object with the list of dataframe added
#' @export
import_exon_annotation <- function(object, annotation_file, event_type){
  annotation_df <- as.data.frame(fread(annotation_file))
  # maek a list for easy access
  rownames(annotation_df) <- annotation_df$exon
  annotation_df$exon <- NULL

  object@exon_annotation[[event_type]] <- annotation_df

  return(object)
}






#' Main import function to import MATS information into object
#'
#' This function reads in the file and add it to the object
#' MATS result is included to perform differential splicing analysis
#' @param object the object
#' @param MATS_file MATS file from processing pipeline
#' @param event_type event_type
#' @return object with the list of dataframe added
#' @export
import_MATS <- function(object, MATS_file, event_type){
  MATS_df <- as.data.frame(fread(MATS_file))
  # maek a list for easy access
  rownames(MATS_df) <- MATS_df$V1
  MATS_df$V1 <- NULL

  object@MATS_list[[event_type]] <- MATS_df

  return(object)
}
