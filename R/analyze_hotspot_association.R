##### SUBFUNCTIONS ######
#' Calculate number of associated AS events for each hotpsot in event type
#'
#'
#' @param input_lm_df input LM df from the pipeline
#' @param event_type SE, A3SS or A5SS
#' @return count df
#' @export
calculate_sig_AS_event_per_HS <- function(input_lm_df, event_type){
  # calculate number of significant AS events per hotspot
  sig_AS_event_per_HS <- as.data.frame(rowSums(input_lm_df > 0))
  colnames(sig_AS_event_per_HS) <- event_type
  sig_AS_event_per_HS$hs_event <- rownames(sig_AS_event_per_HS)
  return(sig_AS_event_per_HS)
}



#' Classify hotspot mutations
#'
#'
#' @param plot_df hotspot plot df
#' @return mutation class df determined for each hotspot
#' @export
classify_hotspot_mutation <- function(plot_df){
  plot_df$hotspot <- rownames(plot_df)
  # get list of mutation classes
  hotspot_list <- unlist(lapply(plot_df$hotspot, function(x) {
    elements <- strsplit(x, ":")[[1]][1:5]
    recombined_string <- sprintf("%s:%s:%s-%s", elements[2], elements[3], elements[4], elements[5])
    return(recombined_string)
  }))
  mutation_class_list_raw <- lapply(plot_df$hotspot, function(x) strsplit(x, ":")[[1]][6])
  mutation_class_list <- unique(unlist(lapply(unlist(lapply(mutation_class_list_raw,
                                                            function(y) strsplit(y, "\n"))),
                                              function(x) strsplit(x, " ")[[1]][1])))

  # generate a dataframe for mutation class count for each hotspot
  mutation_class_heatmap_df <- data.frame(matrix(0,
                                                 nrow = length(hotspot_list),
                                                 ncol = length(mutation_class_list)))
  rownames(mutation_class_heatmap_df) <- hotspot_list
  colnames(mutation_class_heatmap_df) <- mutation_class_list

  # add elements
  dummy <- lapply(seq(length(rownames(mutation_class_heatmap_df))), function(index){
    # print(rownames(mutation_class_heatmap_df)[index])
    x <- mutation_class_list_raw[[index]][1]
    elements <- strsplit(x, "\n")[[1]]
    lapply(seq(length(elements)), function(y_index){


      y <- elements[y_index]
      y_element <- strsplit(y, "\\s+")[[1]]
      mut_class <- y_element[1]
      mut_count <- y_element[2]
      mutation_class_heatmap_df[index, mut_class] <<- as.numeric(mut_count)
      # print(y)
      return()
    })
  })
  return(mutation_class_heatmap_df)
}



#' Subset the summary
#'
#'
#' @param object SPARKS object
#' @param genes_of_interest Genes that you want subset for in terms of hotspot
#' @return Subsetted summary dataframe for downstream functions
#' @export
subset_genes_of_interest_from_summary <- function(object, genes_of_interest){
  # genes of interest are to select few examples that we may be interested in
  if (length(genes_of_interest) == 0){  # if no ones selected
    genes_of_interest <- unique(extract_gene_symbols(rownames(object@summary_asso_df)))
  }

  # check if genes of interest is not empty
  if (length(genes_of_interest) == 0){  # if no gene is found: then it means it is empty and gonna break downstream
    stop("No events for hot key plot found - please revise gene list")
  }

  # select summary_asso_df rows that are in genes of interest from
  summary_df <- object@summary_asso_df
  summary_df_select <- summary_df[extract_gene_symbols(rownames(summary_df)) %in% genes_of_interest, ]

  # check if sumamry df is not empty
  if (dim(summary_df_select)[1] == 0){  # if no gene is found: then it means it is empty and gonna break downstream
    stop("No summary info for hot key plot found - please revise gene list")
  }
  return(summary_df_select)
}



#' Query list of AS events associated for given Mutation event
#'
#' This function generates heatmap plot - x = hotspot mutation status, y = PSI
#'
#' @param object SPARKS object
#' @param mutation_event Mutation event
#' @return list of AS events associated for given Mutation event
#' @export
query_associated_AS_events <- function(object, mutation_event){
  event_types <- c("SE", "A3SS", "A5SS")

  sig_event_df <- do.call(rbind, lapply(event_types, function(event_type){
    # extract data
    event_pval <- object@linear_model_pval[[event_type]][mutation_event, ]
    event_beta <- object@linear_model_beta[[event_type]][mutation_event, ]

    # combine the data
    combined_info <- as.data.frame(t(rbind(event_pval, event_beta)))
    colnames(combined_info) <- c("pval", "beta")

    # add event type for downstream functions
    combined_info$event_type <- event_type

    # keep only relevant information
    filtered_info <- combined_info[combined_info$pval > -log10(0.05), ]

    # drop NA if the rowname is missing
    filtered_info <- filtered_info[rowSums(is.na(filtered_info)) == 0, ]

    return(filtered_info)
  }))
  return(sig_event_df)
}





#' Plot function for hot key plot
#'
#' This function generates hot key plot - row = hotspot, col = mutation class (e.g. Missense mutation, Nonsense mutation)
#'
#' @param mutation_class_heatmap_df_sorted Mutation class heatmap data
#' @return Hotkey plot
#' @export
generate_hotkey_plot_for_hotspot_mutation_type <- function(mutation_class_heatmap_df_sorted){

  # generate a heatmap with mutation class info for annotation
  mutation_class_heatmap_df_sorted$hotspot <- rownames(mutation_class_heatmap_df_sorted)
  colnames(mutation_class_heatmap_df_sorted) <- gsub(x = colnames(mutation_class_heatmap_df_sorted), "_", " ")

  mutation_class_heatmap_melt <- melt(mutation_class_heatmap_df_sorted, id.vars = "hotspot")
  colnames(mutation_class_heatmap_melt) <- c("hotspot", "mutation_type", "count")


  # sort the hotspot into factor
  mutation_class_heatmap_melt$hotspot <- factor(mutation_class_heatmap_melt$hotspot,
                                                levels = rev(rownames(mutation_class_heatmap_df_sorted)))
  mutation_class_heatmap_melt$count <- as.numeric(mutation_class_heatmap_melt$count)
  # mutation_class_heatmap_melt$variable <- gsub(x = mutation_class_heatmap_melt$variable, "_", " ")


  p_mut_class_heatmap <- ggplot(mutation_class_heatmap_melt, aes(x = mutation_type,
                                                                 y = hotspot,
                                                                 fill = mutation_type,
                                                                 alpha = log(count))) +
    geom_raster() +
    geom_text(aes(label = count), size = 3) +
    scale_x_discrete(labels = str_wrap(colnames(mutation_class_heatmap_df_sorted), width = 20)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust = 1,
                                     size = 8),
          axis.text.y = element_text(size = 8)) +
    labs(x = "Mutation\nClassification",
         y = "Somatic Mutation Hotspot")
  return(p_mut_class_heatmap)
}



#' Generate table plot with number of associated AS counts per hotspot
#'
#' This function generates heatmap plot - row = hotspot, col = splicing type, element = # AS events associated
#'
#' @param rbp_hs_sum_df_sorted Hotspot assocation summary heatmap data
#' @return Heatmap
#' @export
generate_table_plot_for_num_asso_AS_per_hotspot <- function(rbp_hs_sum_df_sorted){
  # prep the data for plot
  rbp_hs_sum_df_dummy <- rbp_hs_sum_df_sorted
  rbp_hs_sum_df_dummy$hotspot <- rownames(rbp_hs_sum_df_dummy)
  rbp_hs_sum_melt <- melt(rbp_hs_sum_df_dummy)

  colnames(rbp_hs_sum_melt) <- c("hotspot", "spl_type", "count")

  # set factors for the ordering
  rbp_hs_sum_melt$hotspot <- factor(rbp_hs_sum_melt$hotspot,
                                    levels = rev(rownames(rbp_hs_sum_df_dummy)))
  rbp_hs_sum_melt$spl_type <- factor(rbp_hs_sum_melt$spl_type,
                                     levels = c("SE", "A3SS", "A5SS"))

  p_rbp_hs_asso_count_per_type <- ggplot(rbp_hs_sum_melt, aes(x = spl_type, y = hotspot, fill = log10(count))) +
    geom_raster() +
    geom_text(label = rbp_hs_sum_melt$count, size = 3) +
    scale_fill_distiller(palette = "Oranges", direction = 1, na.value = "white") +
    scale_x_discrete(expand = c(0, 0)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")+
    labs(x = "Splicing Type")
  return(p_rbp_hs_asso_count_per_type)
}


#' Generate boxplot for PSI depending on hotspot mutation status
#'
#' This function generates heatmap plot - x = hotspot mutation status, y = PSI
#'
#' @param object SPARKS object
#' @param hotspot_event Mutation event
#' @param AS_event AS event
#' @return boxplot with sample grouped by mutation status
#' @export
generate_PSI_boxplot_for_mutation_group <- function(object, hotspot_event, AS_event){
  # query sample-wise mutation status for the hotspot mutation
  hotspot_values <- object@hotspot_mutation_matrix[hotspot_event, ]

  # query PSI values for the splicing event
  psi_values <- query_PSI_value_for_AS_event(object, AS_event)

  # bind mutation status and PSI value for plotting
  plot_df <- t(rbind(hotspot_values, psi_values[colnames(hotspot_values)]))
  colnames(plot_df) <- c("hotspot", "psi")

  # generate plot
  p <- ggplot(as.data.frame(plot_df), aes(x = factor(hotspot), y = psi, fill = factor(hotspot))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.7) +
    scale_x_discrete(labels = c("0" = "Without hotspot mutation",
                                "1" = "With hotspot mutation")) +
    labs(x = sprintf("%s samples", object@study),
         y = "Percent Spliced in (PSI)") +
    theme(legend.position = "none")
  # TODO - add set color scheme for samples with/without mutations
  return(p)
}



#' Generate boxplot for EExpression depending on hotspot mutation status
#'
#' This function generates heatmap plot - x = hotspot mutation status, y = Expression
#'
#' @param object SPARKS object
#' @param hotspot_event Mutation event
#' @param AS_event AS event - the gene would be selected by gene symbol
#' @return boxplot with sample grouped by mutation status
#' @export
generate_EXP_boxplot_for_mutation_group  <- function(object, hotspot_event, AS_event,
                                                     exp_df = exp_df_clean,
                                                     hotspot_status = hotspot_status_df){
  # extract mutation status
  hotspot_values <- object@hotspot_mutation_matrix[hotspot_event, ]

  # query expression value
  exp_values <- extract_exp_value_for_AS_event(object, AS_event)

  # bind mutation status and PSI value for plotting
  plot_df <- t(rbind(hotspot_values[colnames(exp_values)], exp_values))
  colnames(plot_df) <- c("hotspot", "exp")

  p <- ggplot(as.data.frame(plot_df), aes(x = factor(hotspot), y = log10(exp), fill = factor(hotspot))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.7) +
    scale_x_discrete(labels = c("0" = "Without hotspot mutation",
                                "1" = "With hotspot mutation")) +
    labs(x = "TCGA-SKCM samples",
         y = "log10(Expression) (TPM)") +
    theme(legend.position = "none")

  return(p)
}



##### MAIN FUNCTIONS #####
#' Generate Summary table for each hotspot assocation
#'
#'
#' @param object SPARKS object
#' @return SPARKS object with the summary df added
#' @export
generate_summary_table_sig_AS_per_HS <- function(object) {
  print("Calculationg summary data for AS association in each hotspot")
  # calculate summary stats
  sig_count_SE <- calculate_sig_AS_event_per_HS(object@linear_model_pval$SE, "SE")
  sig_count_A3SS <- calculate_sig_AS_event_per_HS(object@linear_model_pval$A3SS, "A3SS")
  sig_count_A5SS <- calculate_sig_AS_event_per_HS(object@linear_model_pval$A5SS, "A5SS")

  # merge stats
  sig_stats_df <- full_join(full_join(sig_count_A3SS,
                                      sig_count_A5SS,
                                      by = "hs_event"),
                            sig_count_SE,
                            by = "hs_event")
  rownames(sig_stats_df) <- sig_stats_df$hs_event
  sig_stats_df$hs_event <- NULL

  # fill NA with 0 - the value count should be 0 if the rownames is not in the df
  sig_stats_df[is.na(sig_stats_df)] <- 0

  # remove events that have 0 significant association for downstream analysis
  sig_stats_df_filtered <- sig_stats_df[rowSums(sig_stats_df) > 0, ]

  object@summary_asso_df <- sig_stats_df_filtered
  print("Summary table added - look for summary_asso_df")
  return(object)
}



#' Wrapper for generating hotkey plot for summary
#'
#' @param object SPARKS object
#' @param genes_of_interest Genes you are interested in (hotspot)
#' @param event_group Group of mutations (e.g. Missense mutation) to limit the plot
#' @return hotkey plot for each mutation
#' @export
generate_hotspot_hotkey_plot <- function(object, genes_of_interest = c(), event_group = c()){
  # subset the summary for only interesting genes
  summary_df_select <- subset_genes_of_interest_from_summary(object, genes_of_interest)

  # generate hotkey df for plotting
  hotkey_df <- classify_hotspot_mutation(summary_df_select)
  if (length(event_group) != 0){ # if the event group is selected
    hotkey_df_select <- hotkey_df[hotkey_df[, colnames(hotkey_df) %in% event_group] > 0, ]
    hotkey_df_select <- hotkey_df_select[, colSums(hotkey_df_select) > 0, drop = FALSE]
  } else {
    hotkey_df_select <- hotkey_df
  }

  # generate hotkey plot
  p_hotkey <- generate_hotkey_plot_for_hotspot_mutation_type(hotkey_df_select)

  return(p_hotkey)

}



#' Wrapper for generating heatmap plot for summary
#'
#' @param object SPARKS object
#' @param genes_of_interest Genes you are interested in (hotspot)
#' @param event_group Group of mutations (e.g. Missense mutation) to limit the plot
#' @return summary heatmap for each mutation and associated AS events
#' @export
generate_hotspot_summary_table_plot <- function(object, genes_of_interest = c(), event_group = c()){
  # subset the summary for only interesting genes
  summary_df_select <- subset_genes_of_interest_from_summary(object, genes_of_interest)

  # generate hotkey df for plotting
  hotkey_df <- classify_hotspot_mutation(summary_df_select)
  if (length(event_group) != 0){ # if the event group is selected
    rownames(hotkey_df) <- rownames(summary_df_select)
    mut_event_select <- hotkey_df[hotkey_df[, colnames(hotkey_df) %in% event_group] > 0, ]
    summary_df_select <- summary_df_select[rownames(mut_event_select), ]
  }

  # generate hotkey plot
  p_summary_table <- generate_table_plot_for_num_asso_AS_per_hotspot(summary_df_select)

  return(p_summary_table)

}











# generate sample table for sashimi plot
generate_sample_table_for_rmats2sashimiplot <- function(object, mutation_event, output_dir = "./"){
  # extract mutation status matrix
  mutation_status <- as.data.frame(t(object@hotspot_mutation_matrix[mutation_event, ]))
  colnames(mutation_status) <- c("mutation_status")

  # extract gene symbol for explicit grouping information
  mutation_gene_symbol <- extract_gene_symbols(mutation_event)

  # add grouping
  mutation_status$group <- ifelse(mutation_status$mutation_status == 1,
                                  sprintf("Samples_with_%s_mutation",
                                          mutation_gene_symbol),
                                  sprintf("Samples_without_%s_mutation",
                                          mutation_gene_symbol))
  # add sample name
  mutation_status$sample <- paste0(output_dir, rownames(mutation_status), ".bam")

  # sort it out for rmats2sashimiplot
  mutation_status_sorted <- mutation_status[, c("sample", "group")]

  return(mutation_status_sorted)
}

