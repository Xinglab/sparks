usethis::use_package("R.utils")
usethis::use_package("data.table")
usethis::use_package("dplyr")
usethis::use_package("reshape2")
usethis::use_package("pbmcapply")
usethis::use_package("stringr")
usethis::use_package("ggplot2")
usethis::use_package("RColorBrewer")
usethis::use_package("fgsea")

# SPARKS class
#' @export
SPARKS <- setClass("SPARKS",
                  slots = c(
                    study = "character",
                    num_sample = "numeric",
                    num_mutation = "numeric",
                    sample_chart = "data.frame",
                    # splicing data
                    psi_df = "list",
                    MATS_list = "list",
                    SPARKS_analysis_result = "list",
                    exon_annotation = "list",  # this annotates whether exon AS is event is known
                    # expression data
                    exp_df = "data.frame"
                  ))


##### SUBFUNCTIONS #####
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





#' Main import function to import expression df into object
#'
#' This function reads in the file and filters the data and add it to the object
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
