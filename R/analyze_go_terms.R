generate_go_plot <- function(GO_result, cutoff, aspect_ratio = 0.2, bar_width = 0.7, x_axis.tilt = 0, zero_filter = -0.01){
  library(ggplot2)
  significant_hits <- subset(GO_result, corrected_score < cutoff & corrected_score > zero_filter)
  significant_hits_two <- significant_hits[order(significant_hits$GO_process),]
  significant_hits_two$order = c(1:length(significant_hits_two$Term))
  ggplot(significant_hits_two, aes(x = reorder(Term, order), y = -log10(corrected_score), fill = GO_process)) +
    geom_col(width = bar_width)+
    # facet_grid(.~GO_process) +
    theme(axis.text.x = element_text(hjust=0)) +
    xlab("GO Terms") +
    ylab("-log10(p value)") +
    theme(aspect.ratio = aspect_ratio)+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_x_discrete(position = "top")+
    theme(panel.background = element_blank())+
    theme(axis.text.x=element_text(angle = x_axis.tilt, hjust=0))+
    coord_flip()
}

save_GO_result <- function(GO_result_table, output_dir, output_prefix){
  output_file = paste0(output_dir, "/", "GO_result_table_", output_prefix, ".txt")
  write.table(GO_result_table, output_file, sep = '\t')
  print(paste("The file is saved at:", output_file))
}

generate_go_plot_subset <- function(GO_result, cutoff, aspect_ratio = 0.2, bar_width = 0.7, x_axis.tilt = 0 ){
  for (process in c("BP", "CC", "MF")){
    library(ggplot2)
    significant_hits <- subset(GO_result, corrected_score < cutoff & corrected_score > 0 & GO_process == process)
    significant_hits_two <- significant_hits[order(significant_hits$GO_process),]
    significant_hits_two$order = c(1:length(significant_hits_two$Term))
    ggplot(significant_hits_two, aes(x = reorder(Term, order), y = -log10(corrected_score), fill = GO_process)) +
      geom_col(width = bar_width)+
      # facet_grid(.~GO_process) +
      theme(axis.text.x = element_text(hjust=0)) +
      xlab("GO Terms") +
      ylab("-log10(p value)") +
      theme(aspect.ratio = aspect_ratio)+
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      scale_x_discrete(position = "top")+
      theme(panel.background = element_blank())+
      theme(axis.text.x=element_text(angle = x_axis.tilt, hjust=0))+
      coord_flip() +
      ggsave(filename = paste0("/Users/harryyang/Downloads/go_Plot_", process, ".png"))
  }

}

generate_go_plot_top <- function(GO_result, aspect_ratio = 0.2, bar_width = 0.7, x_axis.tilt = 0, num_top = 10){
  total_sig_hits <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(total_sig_hits) <- colnames(GO_result)
  for (process in c("BP", "CC", "MF")){
    library(ggplot2)
    significant_hits <- subset(GO_result, corrected_score > 0 & GO_process == process )
    significant_hits_top <- significant_hits[1:num_top,]
    print(significant_hits_top)
    total_sig_hits <- rbind(total_sig_hits, significant_hits_top)
  }
  print(total_sig_hits)
  significant_hits_two <- total_sig_hits[order(total_sig_hits$GO_process),]
  significant_hits_two$order = c(1:length(significant_hits_two$Term))
  ggplot(significant_hits_two, aes(x = reorder(Term, order), y = -log10(corrected_score), fill = GO_process)) +
    geom_col(width = bar_width)+
    # facet_grid(.~GO_process) +
    theme(axis.text.x = element_text(hjust=0)) +
    xlab("GO Terms") +
    ylab("-log10(p value)") +
    theme(aspect.ratio = aspect_ratio)+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_x_discrete(position = "top")+
    theme(panel.background = element_blank())+
    theme(axis.text.x=element_text(angle = x_axis.tilt, hjust=0))+
    coord_flip() +
    ggsave(filename = paste0("/Users/harryyang/Downloads/go_Plot_", process, ".png"))


}

generate_go_plot_dev <- function(GO_result, cutoff = 0.05, aspect_ratio = 1.8, bar_width = 0.6, x_axis.tilt = 0, num_top = 10,
                                 text.size = 8){
  total_sig_events <- data.frame(matrix(ncol = 8, nrow = 0)) # generate sig hits data frame for top events.
  colnames(total_sig_events) <- colnames(GO_result)

  for (process in c("BP", "CC", "MF")){
    library(ggplot2)
    significant_hits <- subset(GO_result, GO_process == process ) # subset events for each process
    # TODO - possibly this is the point where we can add the sort
    significant_hits <- significant_hits[order(significant_hits$corrected_score),]
    significant_hits_top <- significant_hits[1:num_top,] # select top number (10?)
    significant_hits_top_after_cutoff <- subset(significant_hits_top, corrected_score < cutoff) # select events below the cutoff p-value
    # print(significant_hits_top)
    total_sig_events <- rbind(total_sig_events, significant_hits_top_after_cutoff)
  }
  print("Top events:")
  print(total_sig_events) # print selected events

  significant_events_ordered <- total_sig_events[order(total_sig_events$GO_process),]
  significant_events_ordered$order = c(1:length(significant_events_ordered$Term)) # put order so that it will be in the order of processes


  # replace GO term code with meaningful numbers
  go_term_new <- unlist(apply(significant_events_ordered, 1, function(go_entry){
    print(go_entry)
    go_term_trimmed <- strsplit(go_entry["Term"], " \\(GO")[[1]][1]
    new_go_term <- sprintf("%s (%i/%i)",
                           go_term_trimmed,
                           as.numeric(go_entry['interest_hit']),
                           as.numeric(go_entry['term_gene_numbers']) )

  }))

  significant_events_ordered$Term_new <- go_term_new


  p <- ggplot(significant_events_ordered, aes(x = reorder(Term_new, -order), y = -log10(corrected_score), fill = GO_process)) +
    geom_col(width = bar_width)+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 50, simplify = FALSE), paste, collapse="\n"), position = "top") +
    theme(axis.text.x = element_text(hjust = 0,
                                     face = "bold",
                                     angle = x_axis.tilt),
          axis.text.y = element_text(size = text.size, vjust = 0.5),
          aspect.ratio = aspect_ratio,
          panel.background = element_blank(),
          legend.position = "bottom")+
    labs(x = "GO Terms", y = "-log10(p value)", fill = "GO Processes") +
    coord_flip()
  return(p)
}

generate_go_plot_top_cutoff <- function(GO_result, cutoff = 0.05, aspect_ratio = 1.8, bar_width = 0.6, x_axis.tilt = 0, num_top = 10){
  total_sig_hits <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(total_sig_hits) <- colnames(GO_result)
  for (process in c("BP", "CC", "MF")){
    library(ggplot2)
    significant_hits <- subset(GO_result, corrected_score > 0 & GO_process == process )
    significant_hits_top <- significant_hits[1:num_top,]
    significant_hits_top_cutoff <- subset(significant_hits_top, corrected_score < cutoff)
    # print(significant_hits_top)
    total_sig_hits <- rbind(total_sig_hits, significant_hits_top_cutoff)
  }
  print(total_sig_hits)
  significant_hits_two <- total_sig_hits[order(total_sig_hits$GO_process),]
  significant_hits_two$order = c(1:length(significant_hits_two$Term))
  ggplot(significant_hits_two, aes(x = reorder(Term, order), y = -log10(corrected_score), fill = GO_process)) +
    geom_col(width = bar_width)+
    # facet_grid(.~GO_process) +
    theme(axis.text.x = element_text(hjust=0)) +
    xlab("GO Terms") +
    ylab("-log10(p value)") +
    theme(aspect.ratio = aspect_ratio)+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_x_discrete(position = "top")+
    theme(panel.background = element_blank())+
    theme(axis.text.x=element_text(angle = x_axis.tilt, hjust=0))+
    coord_flip() +
    ggsave(filename = paste0("/Users/harryyang/Downloads/go_Plot_", process, ".png"))


}


get_enrichr_GO <- function(all_gene_list, significant_genes){
  GO_processes = c("BP", "CC", "MF")
  merged_GO_result <- list()
  for (process in GO_processes){
    recalculated_GO_result = recalculate_GO(all_gene_list, significant_genes, process)
    merged_GO_result <- rbind(recalculated_GO_result, merged_GO_result)
  }
  merged_GO_result

}
recalculate_GO <- function(all_gene_list, significant_genes, GO_process = "BP") {
  library(enrichR)
  if (GO_process == "BP") {
    go_db = "GO_Biological_Process_2021"
  }  else if (GO_process == "CC"){
    go_db = "GO_Cellular_Component_2021"
  }  else if (GO_process == "MF"){
    go_db = "GO_Molecular_Function_2021"
  }
  all_gene_GO <- enrichr(as.vector(all_gene_list), databases = go_db)
  sig_gene_GO <- enrichr(as.vector(significant_genes), databases = go_db)

  bg_numbers <- as.data.frame(cbind(all_gene_GO[[go_db]]$Term,
                                    sapply(strsplit(as.character(all_gene_GO[[go_db]]$Overlap), '/'), '[', 1),
                                    sapply(strsplit(as.character(all_gene_GO[[go_db]]$Overlap), '/'), '[', 2)),
                              stringsAsFactors = FALSE)

  colnames(bg_numbers) <- c("Term", "Overlap_number", "Term_gene_number")
  rownames(bg_numbers) <- bg_numbers$Term

  sig_numbers <- as.data.frame(cbind(sig_gene_GO[[go_db]]$Term,
                                     sapply(strsplit(as.character(sig_gene_GO[[go_db]]$Overlap), '/'), '[', 1),
                                     sapply(strsplit(as.character(sig_gene_GO[[go_db]]$Overlap), '/'), '[', 2)), stringsAsFactors = FALSE)
  colnames(sig_numbers) <- c("Term", "Overlap_number", "Term_gene_number")
  rownames(sig_numbers) <- sig_numbers$Term

  corrected_score <- sapply(sig_numbers$Term, function(x) get_hypergeometric_values(x, sig_numbers, bg_numbers, all_gene_list, significant_genes))
  names(corrected_score) <- sig_numbers$Term

  # remove NA
  corrected_score <- corrected_score[!is.na(corrected_score)]

  # FDR Correction
  corrected_score_after_fdr <- p.adjust(corrected_score, method = "fdr")

  # formatting
  corrected_numbers <- cbind(as.numeric(as.character(sig_numbers[names(corrected_score_after_fdr),]$Overlap_number)),
                             as.numeric(as.character(bg_numbers[names(corrected_score_after_fdr),]$Overlap_number)),
                             as.numeric(as.character(bg_numbers[names(corrected_score_after_fdr),]$Term_gene_number)))
  colnames(corrected_numbers) <- c("Sig_event_num", "Background_event_num", "Term_gene_number")
  rownames(corrected_numbers) <- sig_numbers[names(corrected_score_after_fdr),]$Term
  # combined_num <- sapply(corrected_numbers, function(x) as.character(x["Sig_event_num"]/x["Background_event_num"]))


  corrected_genes <- sig_gene_GO[[go_db]][sig_gene_GO[[go_db]]$Term %in% rownames(corrected_numbers),]$Genes

  info_table <- cbind(rownames(corrected_numbers),
                      corrected_numbers,
                      data.frame(corrected_score),
                      data.frame(corrected_score_after_fdr),
                      corrected_genes)
  # print(head(info_table))
  colnames(info_table) <- c("Term", "interest_hit", "background_hit", "term_gene_numbers",
                            "original_score", "corrected_score", "Genes")
  info_table$GO_process <- GO_process
  info_table
}

get_hypergeometric_values <- function (term, sig_num, bg_num, all_gene_list, sig_gene_list)  {
  term_bg <- as.numeric(as.character(bg_num[term,]$Overlap_number))
  # print(term_bg)
  term_bg[is.na(term_bg)] <- 0

  if (term_bg >= 5){
    term_obs <- as.numeric(as.character(sig_num[term,]$Overlap_number))
    term_bg <- as.numeric(as.character(bg_num[term,]$Overlap_number))
    total_bg <- length(all_gene_list)
    total_obs <- length(sig_gene_list)
    # pval=hypergeom.sf(term_obs,tot_bg,term_bg,tot_obs)
    # sf(x, M, n, N, loc=0)
    updated_p_value = 1 - phyper(term_obs, term_bg, total_bg - term_bg, total_obs, lower.tail = TRUE)
    updated_p_value
  }
  else{
    # print("WRONG ")
    NA
  }
}


##### MAIN FUNCTION #####
generate_GO_plot <- function(object, mutation_event){
  gene_list <- lapply(names(object@linear_model_pval), function(spl_event_type){

    mutation_asso_event_df <- object@linear_model_pval[[spl_event_type]][mutation_event, ]
    mutation_asso_event <- colnames(mutation_asso_event_df)[mutation_asso_event_df > 0]
    mutation_asso_gene <- unique(extract_gene_symbols(mutation_asso_event))

  })
  unique_asso_gene <- unique(unlist(gene_list))

  # extract background gene list
  background_gene_list <- lapply(names(object@linear_model_pval), function(spl_event_type){

    all_AS <- colnames(object@linear_model_pval[[spl_event_type]])
    all_gene <- unique(extract_gene_symbols(all_AS))

  })
  background_gene <- unique(unlist(background_gene_list))
  print(unique_asso_gene)
  # perform GO analysis
  GO_table_super_AS <- get_enrichr_GO(background_gene, unique_asso_gene)
  p_GO <- generate_go_plot_dev(GO_table_super_AS, cutoff = 0.05)
  return(p_GO)
}
