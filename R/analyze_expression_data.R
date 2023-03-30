
#' Query expression data from SPARKS object
#'
#' @param object SPARKS object
#' @param gene_list list of genes you are interested in
#'
#' @return Dataframe of expression values for the select genes
#' @export
query_expression_data_for_gene_set <- function(object,
                                               gene_list){
  expression_df <- object@exp_df

  trimmed_gene_name <- extract_gene_symbols(rownames(expression_df))

  target_exp_df <- expression_df[trimmed_gene_name %in% gene_list, ]
  # add the gene for melting for downstream plotting
  target_exp_df$gene <- extract_gene_symbols(rownames(target_exp_df))
  rownames(target_exp_df) <- NULL
  return(target_exp_df)
}





