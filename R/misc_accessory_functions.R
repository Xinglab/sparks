#' Query PSI value for AS event
#'
#'
#' @param object SPARKS object
#' @param event_type AS event
#' @return psi values for the AS event
#' @export
query_PSI_value_for_AS_event <- function(object, AS_event) {
  # figure out event type and correct PSI df
  event_type <- strsplit(AS_event, ":")[[1]][length(strsplit(AS_event, ":")[[1]])]
  if (event_type == "A3SS"){
    psi_df <- object@psi_df$A3SS
  } else if (event_type == "A5SS"){
    psi_df <- object@psi_df$A5SS
  } else if (event_type == "SE"){
    psi_df <- object@psi_df$SE
  }
  psi_values <- psi_df[AS_event, ]
  return(psi_values)
}




#' Query EXP value for AS event
#'
#'
#' @param object SPARKS object
#' @param event_type AS event
#' @return EXP values for the AS event
#' @export
extract_exp_value_for_AS_event <- function(object, AS_event) {
  # figure out event type and correct PSI df
  gene_id <- strsplit(strsplit(AS_event, ":")[[1]][1], "\\.")[[1]][1]
  exp_select <- object@exp_df[grep(gene_id, rownames(object@exp_df)),]
  return(exp_select)
}



#' Extract gene symbols from complex annotation
#'
#' @param hotspot_list AS event/hotspot list
#' @return Gene symbol, which is generally in the second place in colon seperated names
#' @export
extract_gene_symbols <- function(hotspot_list){
  # extract genes from hotspots or alt splicing events
  # this functino works for both since the gene names are the second term
  # when colon seperated.
  gene_symbol_list <- unlist(lapply(hotspot_list, function(x) strsplit(x, ":")[[1]][2]))
  return(gene_symbol_list)
}




