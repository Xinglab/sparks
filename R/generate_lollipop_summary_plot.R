#' Generate Lollipop Plot with the Strip for specific Splicing Type
#'
#' @param test_result_list SPARKS analysis result list for SE, A3SS, A5SS
#' @param rbp_of_interest RBPs that would be annotated on the lollipop
#' @param spl_type Splicing Type - SE, A5SS, A3SS
#' @return Lollipop plot with strip on the bottom for tracks
#' @export
generate_strip_lollipop_plot <- function(test_result_list,
                                         rbp_of_interest,
                                         spl_type,
                                         max_score_threshold = 0.8,
                                         color_offset = 0,
                                         score_legend = FALSE,
                                         color = "Set1",
                                         manual_colors = NULL){
  # check validity of input splicing type
  spl_types <- c("SE", "A5SS", "A3SS")
  if (!(spl_type) %in% spl_types){
    stop("Splicing Type is undefined - please check the splice type")
  }

  if (is.null(manual_colors)){
    # manually set colors for RBPs in the rbp_of_interest
    rbp_colors <- RColorBrewer::brewer.pal(9, "Set1")
    if (color_offset > 0){  # shift the elements a bit for better color based on user preference
      rbp_colors <- c(rbp_colors[-seq(color_offset)],
                      rbp_colors[seq(color_offset)])
    }
  } else {  # if colors are manually given
    # validate the input
    if (length(manual_colors) < length(rbp_of_interest)){  # check if the color list is shorter than supplied
      # longer one is allowed so that user can supply an entire set
      stop("Manual colors supplied are shorter than the list of RBPs - please check the input")
    }
    rbp_colors <- manual_colors
  }
  # set names so downstream plots can fit in the color
  # - this depends on the default behavior that the elements past the supplied RBP list would be just NA
  names(rbp_colors) <- rbp_of_interest


  # get test result
  test_result_df <- get_test_results(test_result_list,
                                     rbp_of_interest,
                                     spl_type,
                                     max_score_threshold)

  ## generate lollipop plot
  lolli_point_df <- subset(test_result_df,
                           interest == "interest")

  # generate new labels for lollipop plot annotation
  new_labels <- apply(lolli_point_df, 1, function(entry){
    # NOTE - flatten the rank
    # - the rank can be in 0.5 increment if there are tie
    # - the tie happens in the menial places so doesn
    label_string <- sprintf("%s: %s - %s %s, ES = %s",
                            floor(as.numeric(entry['plot_rank'])),
                            entry['RBP'],
                            entry['cell_line'],
                            entry['strategy'],
                            round(as.numeric(entry['plot_scores']), 3))


    # print(entry)
    # label_string <- sprintf("%s - %s %s\nRank: %s\nScore: %s", entry['RBP'], entry['cell_line'], entry['strategy'],  entry['plot_rank'], round(as.numeric(entry['plot_scores']), 3))
    # label_string <- sprintf("%s - %s %s", entry['RBP'], entry['cell_line'], entry['strategy'])
    # label_string <- sprintf("%s: %s - %s %s\nScore: %s", entry['plot_rank'], entry['RBP'], entry['cell_line'], entry['strategy'],   round(as.numeric(entry['plot_scores']), 2))
    # label_string <- sprintf("%s: %s - %s %s", entry['plot_rank'], entry['RBP'], entry['cell_line'], entry['strategy'])
    return(label_string)
  })
  lolli_point_df$new_label <- new_labels

  # generate the plot
  p_lollipop <-
    ggplot(lolli_point_df,
           aes(x = plot_rank)) +
    geom_segment(aes(yend = 1,  # generate the black line for lollipop
                     y = 0,
                     xend = plot_rank,
                     x = plot_rank),
                 alpha = 1) +
    geom_point(aes(y = 1,  # generate the lollipop
                   color = RBP,
                   fill = flat_score),
               size = 5, shape = 21) +
    ggrepel::geom_text_repel(aes(label = new_labels,  # add the labels to the lollipop
                                 y = 1,
                                 color = RBP),
                             force_pull   = 0, # do not pull toward data points
                             nudge_y      = 0.1,
                             direction    = "x",
                             angle        = 90,
                             hjust        = 0,
                             segment.size = 0.4,
                             max.iter = 1e4,
                             max.time = 1,
                             # fontface = 'bold',
                             family = "Arial",
                             size = 3) +
    scale_x_continuous(limits = c(0, max(test_result_df$plot_rank)),
                       expand = expansion(add = c(10, 10)),
                       breaks = c(seq(0,  # breaks at the beginning
                                      floor(max(test_result_df$plot_rank) / 200)) * 200,  # breaks at every 200
                                  max(test_result_df$plot_rank))) +  # breaks at the end
    scale_y_continuous(breaks = NULL,
                       expand = expansion(mult = c(0, 0.5))) +  # make space for the annotation
    coord_cartesian(ylim = c(0.75, 1.5)) +  # crop the bottom part that was necessary for the lollipop stick
    scale_color_manual(values = rbp_colors,  # set colors for the RBPs in lollipop and labels
                       guide = "none") +
    scale_fill_distiller(palette = "PiYG",  # manually set fill for enrichment score
                         limits = c(-max_score_threshold, max_score_threshold),
                         na.value = "black",
                         breaks = c(0.8, 0.4, 0, -0.4, -0.8),
                         labels = c(">0.8", 0.4, 0, -0.4, "<-0.8"),
                         values = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(fill = "Enrichment Score",
         x = "Enrichment Score Rank") +
    theme(axis.title.y = element_blank(),  # remove axis titles because it is unnecessary
          axis.title.x = element_blank()) +
    theme(panel.background = element_blank(),  # clean out the plot area
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black",  # generate border for better presentation
                                      fill = NA,
                                      size = 1))

  # generate strip plot
  p_strip <-
    ggplot(test_result_df, aes(x = plot_rank, y = spl_type)) +
    geom_tile(aes(fill = flat_score),
              width = 1.1) +  # this width may be necessary to remove weird white lines
    scale_x_continuous(limits = c(0, max(test_result_df$plot_rank)),  # manually set x axis scales to match the strip plot
                       expand = expansion(add = c(10, 10)),
                       breaks = c(seq(0,  # breaks at the begining
                                      floor(max(test_result_df$plot_rank) / 150)) * 150,  # breaks at every 150 experiments
                                  max(test_result_df$plot_rank))) +  # breaks at the end
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_distiller(palette = "PiYG",  # manually set fill for enrichment score
                         limits = c(-max_score_threshold, max_score_threshold),
                         na.value = "black",
                         breaks = c(0.8, 0.4, 0, -0.4, -0.8),
                         labels = c(">0.8", 0.4, 0, -0.4, "<-0.8"),
                         values = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(fill = "Enrichment Score",
         x = "Enrichment Score Rank") +
    guides(fill = guide_colorbar(reverse = T)) +  # reverse so that it would match the strip plot
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),  # change legend text size
          legend.title = element_text(size = 10)) +
    theme(axis.title.y = element_blank(), # change axis title stuffs
          axis.title.x = element_text(size = 10)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1))

  if (!score_legend){  # remove legend for stacked plots
    p_strip <- p_strip + theme(legend.position = "none")
    # combine the plots for final
    p_combined <-
      cowplot::plot_grid(p_lollipop +  # remove axis and legends because it will be annotated in the strip plot
                           theme(legend.position = "none") +
                           theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title = element_blank()),
                         p_strip,
                         ncol = 1,
                         align = "v",
                         axis = "lr",
                         rel_heights = c(5, 1))  #note the ratio is a bit different
  } else {  # remove legend for stacked plots
    # combine the plots for final
    p_combined <-
      cowplot::plot_grid(p_lollipop +  # remove axis and legends because it will be annotated in the strip plot
                           theme(legend.position = "none") +
                           theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title = element_blank()),
                         p_strip,
                         ncol = 1,
                         align = "v",
                         axis = "lr",
                         rel_heights = c(2, 1))
  }


  return(p_combined)
}

#' Generate Lollipop Plot with the Strip for specific Splicing Type
#'
#' @param test_result_list SPARKS analysis result list for SE, A3SS, A5SS
#' @param rbp_of_interest RBPs that would be annotated on the lollipop
#' @param spl_type Splicing Type - SE, A5SS, A3SS
#' @return Lollipop plot with strip on the bottom for tracks
#' @export
generate_strip_lollipop_plot_dev <- function(test_result_list,
                                             rbp_of_interest,
                                             spl_type,
                                             max_score_threshold = 0.8,
                                             color_offset = 0,
                                             # score_legend = TRUE,
                                             color = "Set1",
                                             manual_colors = NULL){
  # check validity of input splicing type
  spl_types <- c("SE", "A5SS", "A3SS")
  if (!(spl_type) %in% spl_types){
    stop("Splicing Type is undefined - please check the splice type")
  }

  if (is.null(manual_colors)){
    # manually set colors for RBPs in the rbp_of_interest
    rbp_colors <- RColorBrewer::brewer.pal(9, "Set1")
    if (color_offset > 0){  # shift the elements a bit for better color based on user preference
      rbp_colors <- c(rbp_colors[-seq(color_offset)], rbp_colors[seq(color_offset)])
    }
  } else {  # if colors are manually given
    # validate the input
    if (length(manual_colors) < length(rbp_of_interest)){  # check if the color list is shorter than supplied
      # longer one is allowed so that user can supply an entire set
      stop("Manual colors supplied are shorter than the list of RBPs - please check the input")
    }
    rbp_colors <- manual_colors
  }
  # set names so downstream plots can fit in the color
  # - this depends on the default behavior that the elements past the supplied RBP list would be just NA
  names(rbp_colors) <- rbp_of_interest


  # get test result
  test_result_df <- get_test_results(test_result_list, rbp_of_interest, spl_type, max_score_threshold)


  ## generate lollipop plot
  lolli_point_df <- subset(test_result_df, interest == "interest")

  # generate new labels for lollipop plot annotation
  new_labels <- apply(lolli_point_df, 1, function(entry){
    # NOTE - flatten the rank
    # - the rank can be in 0.5 increment if there are tie
    # - the tie happens in the menial places so doesn
    label_string <- sprintf("%s: %s - %s %s, ES = %s", floor(as.numeric(entry['plot_rank'])), entry['RBP'], entry['cell_line'], entry['strategy'],   round(as.numeric(entry['plot_scores']), 3))


    # print(entry)
    # label_string <- sprintf("%s - %s %s\nRank: %s\nScore: %s", entry['RBP'], entry['cell_line'], entry['strategy'],  entry['plot_rank'], round(as.numeric(entry['plot_scores']), 3))
    # label_string <- sprintf("%s - %s %s", entry['RBP'], entry['cell_line'], entry['strategy'])
    # label_string <- sprintf("%s: %s - %s %s\nScore: %s", entry['plot_rank'], entry['RBP'], entry['cell_line'], entry['strategy'],   round(as.numeric(entry['plot_scores']), 2))
    # label_string <- sprintf("%s: %s - %s %s", entry['plot_rank'], entry['RBP'], entry['cell_line'], entry['strategy'])
    return(label_string)
  })
  lolli_point_df$new_label <- new_labels

  # generate the plot
  p_lollipop <-
    ggplot(lolli_point_df, aes(x = plot_rank)) +
    geom_segment(aes(yend = 1,  # generate the black line for lollipop
                     y = 0,
                     xend = plot_rank,
                     x = plot_rank),
                 alpha = 1) +
    geom_point(aes(y = 1,  # generate the lollipop
                   color = RBP,
                   fill = flat_score),
               size = 5, shape = 21) +
    scale_fill_distiller(palette = "PiYG",  # manually set fill for enrichment score
                         limits = c(-max_score_threshold, max_score_threshold),
                         na.value = "black",
                         breaks = c(0.8, 0.4, 0, -0.4, -0.8),
                         labels = c(">0.8", 0.4, 0, -0.4, "<-0.8"),
                         values = c(0, 0.25, 0.5, 0.75, 1)) +
    ggrepel::geom_text_repel(aes(label = new_labels,  # add the labels to the lollipop
                                 y = 1,
                                 color = RBP),
                             force_pull   = 0, # do not pull toward data points
                             nudge_y      = 0.1,
                             direction    = "x",
                             angle        = 90,
                             hjust        = 0,
                             segment.size = 0.4,
                             max.iter = 1e4,
                             max.time = 1,
                             # fontface = 'bold',
                             family = "Arial",
                             size = 3) +
    scale_color_manual(values = rbp_colors,  # set colors for the RBPs in lollipop and labels
                       guide = "none") +
    scale_y_continuous(breaks = NULL,
                       expand = expansion(mult = c(0, 0.5))) +  # make space for the annotation
    coord_cartesian(ylim = c(0.75, 1.5)) +  # crop the bottom part that was necessary for the lollipop stick
    scale_x_continuous(limits = c(0, max(test_result_df$plot_rank)),
                       expand = expansion(add = c(10, 10)),
                       breaks = c(seq(0,  # breaks at the beginning
                                      floor(max(test_result_df$plot_rank) / 200)) * 200,  # breaks at every 200
                                  max(test_result_df$plot_rank))) +  # breaks at the end
    guides(fill = guide_colorbar(reverse = T,
                                 title.position = "top")) +
    labs(fill = "Enrichment Score",
         x = "Enrichment Score Rank") +
    # theme(legend.text = element_text(size = 8),  # change legend text size
    #       legend.title = element_text(size = 10)) +
    theme(axis.title.y = element_blank(),  # remove axis titles because it is unnecessary
          axis.title.x = element_blank()) +
    theme(panel.background = element_blank(),  # clean out the plot area
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black",  # generate border for better presentation
                                      fill = NA,
                                      size = 1)) +
    theme(legend.position = c(0.75, 0.9),
          legend.background = element_rect(fill = alpha("white", 0.5),
                                           # color = "black"
                                           ),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.direction = "horizontal",
          legend.title.align = 0.5,
          legend.margin = margin(t = 5, b = 5, r = 15, l = 15),
          legend.spacing.x = unit(10, "pt"),
          legend.key.width = unit(20, 'pt'),
          legend.key.height = unit(5, "pt"))


  # generate strip plot
  p_strip <-
    ggplot(test_result_df, aes(x = plot_rank, y = spl_type)) +
    geom_tile(aes(fill = flat_score),
              width = 1.1) +  # this width may be necessary to remove weird white lines
    scale_fill_distiller(palette = "PiYG",  # manually set fill for enrichment score
                         limits = c(-max_score_threshold, max_score_threshold),
                         na.value = "black",
                         breaks = c(0.8, 0.4, 0, -0.4, -0.8),
                         labels = c(">0.8", 0.4, 0, -0.4, "<-0.8"),
                         values = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(limits = c(0, max(test_result_df$plot_rank)),  # manually set x axis scales to match the strip plot
                       expand = expansion(add = c(10, 10)),
                       breaks = c(seq(0,  # breaks at the begining
                                      floor(max(test_result_df$plot_rank) / 150)) * 150,  # breaks at every 150 experiments
                                  max(test_result_df$plot_rank))) +  # breaks at the end
    scale_y_discrete(expand = c(0, 0)) +
    labs(fill = "Enrichment Score",
         x = "Enrichment Score Rank") +
    theme(legend.position = "bottom") +
    guides(fill = guide_colorbar(reverse = T)) +  # reverse so that it would match the strip plot
    theme(legend.text = element_text(size = 8),  # change legend text size
          legend.title = element_text(size = 10)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_text(size = 10)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1))


    p_strip <- p_strip +
      theme(legend.position = "none")
    # combine the plots for final
    p_combined <-
      cowplot::plot_grid(p_lollipop +  # remove axis and legends because it will be annotated in the strip plot
                           theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title = element_blank()) +
                           theme(legend.position = "none"),
                         p_strip,
                         ncol = 1,
                         align = "v",
                         axis = "lr",
                         rel_heights = c(5, 1))  #note the ratio is a bit different



  return(p_combined)
}


generate_strip_legend <- function(test_result_list,
                                  rbp_of_interest,
                                  spl_type,
                                  max_score_threshold = 0.8,
                                  color_offset = 0,
                                  score_legend = TRUE){
  # check validity of input splicing type
  spl_types <- c("SE", "A5SS", "A3SS")
  if (!(spl_type) %in% spl_types){
    stop("Splicing Type is undefined - please check the splice type")
  }

  # manually set colors for RBPs in the rbp_of_interest
  rbp_colors <- RColorBrewer::brewer.pal(9, "Set1")
  if (color_offset > 0){  # shift the elements a bit for better color based on user preference
    rbp_colors <- c(rbp_colors[-seq(color_offset)], rbp_colors[seq(color_offset)])
  }
  names(rbp_colors) <- rbp_of_interest

  # get test result
  test_result_df <- get_test_results(test_result_list, rbp_of_interest, spl_type, max_score_threshold)


  # generate strip plot
  p_strip <-
    ggplot(test_result_df, aes(x = plot_rank, y = spl_type)) +
    geom_tile(aes(fill = flat_score),
              width = 1.1) +  # this width may be necessary to remove weird white lines
    scale_fill_distiller(palette = "PiYG",  # manually set fill for enrichment score
                         limits = c(-max_score_threshold, max_score_threshold),
                         na.value = "black",
                         breaks = c(0.8, 0.4, 0, -0.4, -0.8),
                         labels = c(">0.8", 0.4, 0, -0.4, "<-0.8"),
                         values = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(limits = c(0, max(test_result_df$plot_rank)),  # manually set x axis scales to match the strip plot
                       expand = expansion(add = c(10, 10)),
                       breaks = c(seq(0,  # breaks at the begining
                                      floor(max(test_result_df$plot_rank) / 150)) * 150,  # breaks at every 150 experiments
                                  max(test_result_df$plot_rank))) +  # breaks at the end
    scale_y_discrete(expand = c(0, 0)) +
    labs(fill = "Enrichment Score",
         x = "Enrichment Score Rank") +
    theme(legend.position = "bottom") +
    guides(fill = guide_colorbar(reverse = T)) +  # reverse so that it would match the strip plot
    theme(legend.text = element_text(size = 8),  # change legend text size
          legend.title = element_text(size = 10)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_text(size = 10)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1))

  return(cowplot::get_legend(p_strip))
  # return(p_combined)
}
##### AUX FUNCTION #####

get_test_results <- function(test_result_list, rbp_of_interest, spl_type, max_score_threshold){


  test_result_df <- test_result_list[[spl_type]]

  # sort by plot score
  test_result_df$plot_scores <- test_result_df$gsea_score
  test_result_df$S1 <- factor(test_result_df$S1, levels = test_result_df[order(-test_result_df$plot_scores), ]$S1)
  test_result_df$plot_rank <- rank(-test_result_df$plot_scores)

  # flatten the score above the threshold
  test_result_df$flat_score <- test_result_df$plot_scores
  test_result_df$flat_score[test_result_df$flat_score > max_score_threshold] <- max_score_threshold
  test_result_df$flat_score[test_result_df$flat_score < -max_score_threshold] <- -max_score_threshold

  # annotate RBP and flag for highlight
  test_result_df$RBP <- unlist(lapply(as.character(test_result_df$S1), function(x) strsplit(x, "_")[[1]][2]))
  test_result_df$interest <- ifelse(test_result_df$RBP %in% c(rbp_of_interest), "interest", NA)
  # test_result_df[!is.na(test_result_df$interest), ]$plot_scores <- NA  # Using NA to for manual setting

  # annotate other stuffs
  test_result_df$cell_line <- unlist(lapply(as.character(test_result_df$S1), function(x) strsplit(x, "_")[[1]][1]))
  test_result_df$strategy <- unlist(lapply(as.character(test_result_df$S1), function(x) strsplit(x, "_")[[1]][3]))

  return(test_result_df)
}

#' Generate Lollipop Plot with the Strip for specific Splicing Type
#'
#' @param test_result_list SPARKS analysis result list for SE, A3SS, A5SS
#' @param rbp_of_interest RBPs that would be annotated on the lollipop
#' @param spl_type Splicing Type - SE, A5SS, A3SS
#' @return Lollipop plot with strip on the left for tracks
#' @export
generate_strip_lollipop_plot_vertical <- function(test_result_list,
                                                  rbp_of_interest,
                                                  spl_type,
                                                  max_score_threshold = 0.8,
                                                  color_offset = 0,
                                                  score_legend = FALSE,
                                                  color = "Set1",
                                                  manual_colors = NULL){
  # check validity of input splicing type
  spl_types <- c("SE", "A5SS", "A3SS")
  if (!(spl_type) %in% spl_types){
    stop("Splicing Type is undefined - please check the splice type")
  }

  if (is.null(manual_colors)){
    # manually set colors for RBPs in the rbp_of_interest
    rbp_colors <- RColorBrewer::brewer.pal(9, "Set1")
    if (color_offset > 0){  # shift the elements a bit for better color based on user preference
      rbp_colors <- c(rbp_colors[-seq(color_offset)],
                      rbp_colors[seq(color_offset)])
    }
  } else {  # if colors are manually given
    # validate the input
    if (length(manual_colors) < length(rbp_of_interest)){  # check if the color list is shorter than supplied
      # longer one is allowed so that user can supply an entire set
      stop("Manual colors supplied are shorter than the list of RBPs - please check the input")
    }
    rbp_colors <- manual_colors
  }
  # set names so downstream plots can fit in the color
  # - this depends on the default behavior that the elements past the supplied RBP list would be just NA
  names(rbp_colors) <- rbp_of_interest


  # get test result
  test_result_df <- get_test_results(test_result_list,
                                     rbp_of_interest,
                                     spl_type,
                                     max_score_threshold)

  ## generate lollipop plot
  lolli_point_df <- subset(test_result_df,
                           interest == "interest")

  # generate new labels for lollipop plot annotation
  new_labels <- apply(lolli_point_df, 1, function(entry){
    # NOTE - flatten the rank
    # - the rank can be in 0.5 increment if there are tie
    # - the tie happens in the menial places so doesn
    # label_string <- sprintf("%s: %s - %s %s\nES = %s",
    #                         floor(as.numeric(entry['plot_rank'])),
    #                         entry['RBP'],
    #                         entry['cell_line'],
    #                         entry['strategy'],
    #                         round(as.numeric(entry['plot_scores']), 3))
    #
    label_string <- sprintf("%s %s %s\nRank = %s, ES = %s",
                            entry['RBP'],
                            entry['cell_line'],
                            entry['strategy'],
                            floor(as.numeric(entry['plot_rank'])),
                            round(as.numeric(entry['plot_scores']), 3))

    # print(entry)
    # label_string <- sprintf("%s - %s %s\nRank: %s\nScore: %s", entry['RBP'], entry['cell_line'], entry['strategy'],  entry['plot_rank'], round(as.numeric(entry['plot_scores']), 3))
    # label_string <- sprintf("%s - %s %s", entry['RBP'], entry['cell_line'], entry['strategy'])
    # label_string <- sprintf("%s: %s - %s %s\nScore: %s", entry['plot_rank'], entry['RBP'], entry['cell_line'], entry['strategy'],   round(as.numeric(entry['plot_scores']), 2))
    # label_string <- sprintf("%s: %s - %s %s", entry['plot_rank'], entry['RBP'], entry['cell_line'], entry['strategy'])
    return(label_string)
  })
  lolli_point_df$new_label <- new_labels

  # generate the plot
  p_lollipop <-
    ggplot(lolli_point_df,
           aes(y = plot_rank)) +
    geom_segment(aes(xend = 1,  # generate the black line for lollipop
                     x = 0,
                     yend = plot_rank,
                     y = plot_rank),
                 alpha = 1) +
    geom_point(aes(x = 1,  # generate the lollipop
                   color = RBP,
                   fill = flat_score),
               size = 5, shape = 21) +
    ggrepel::geom_text_repel(aes(label = new_labels,  # add the labels to the lollipop
                                 x = 1,
                                 color = RBP),
                             force_pull   = 0, # do not pull toward data points
                             nudge_x      = 0.1,
                             direction    = "y",
                             angle        = 0,
                             hjust        = 0,
                             segment.size = 0.4,
                             max.iter = 1e4,
                             max.time = 1,
                             # fontface = 'bold',
                             family = "Arial",
                             size = 3) +
    scale_y_reverse(limits = rev(c(0, max(test_result_df$plot_rank))),
                    expand = expansion(add = c(10, 10)),
                    breaks = c(1, # breaks at the beginning
                               seq(floor(max(test_result_df$plot_rank) / 200)) * 200,  # breaks at every 200
                               max(test_result_df$plot_rank))) +  # breaks at the end
    scale_x_continuous(breaks = NULL,
                       expand = expansion(mult = c(0, 0.5))) +  # make space for the annotation
    coord_cartesian(xlim = c(0.9, 1.5)) +  # crop the bottom part that was necessary for the lollipop stick
    scale_color_manual(values = rbp_colors[!is.na(names(rbp_colors))]) +  # remove RBP colors with no entries
    scale_fill_distiller(palette = "PiYG",  # manually set fill for enrichment score
                         limits = c(-max_score_threshold, max_score_threshold),
                         na.value = "black",
                         breaks = c(0.8, 0.4, 0, -0.4, -0.8),
                         labels = c(">0.8", 0.4, 0, -0.4, "<-0.8"),
                         values = c(0, 0.25, 0.5, 0.75, 1),
                         guide = "none") +
    labs(fill = "Enrichment Score",
         x = "Enrichment Score Rank") +
    theme(axis.title.y = element_blank(),  # remove axis titles because it is unnecessary
          axis.title.x = element_blank()) +
    theme(panel.background = element_blank(),  # clean out the plot area
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black",  # generate border for better presentation
                                      fill = NA,
                                      size = 1))#+
  # theme(legend.position = c(0.5, 0.15),
  #       legend.background = element_rect(fill = "white", color = "black"),
  #       legend.text = element_text(size = 8),
  #       legend.title = element_text(size = 8)) ### PLOT SIZE - w 700 x h 400

  # generate strip plot
  p_strip <-
    ggplot(test_result_df, aes(y = plot_rank, x = spl_type)) +
    geom_tile(aes(fill = flat_score),
              width = 1.1) +  # this width may be necessary to remove weird white lines
    scale_y_reverse(limits = rev(c(0, max(test_result_df$plot_rank))),  # manually set x axis scales to match the strip plot
                    expand = expansion(add = c(10, 10)),
                    breaks = c(1, # breaks at the begining
                               seq(floor(max(test_result_df$plot_rank) / 150)) * 150,  # breaks at every 150 experiments
                               max(test_result_df$plot_rank))) +  # breaks at the end
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_distiller(palette = "PiYG",  # manually set fill for enrichment score
                         limits = c(-max_score_threshold, max_score_threshold),
                         na.value = "black",
                         breaks = c(0.8, 0.4, 0, -0.4, -0.8),
                         labels = c(">0.8", 0.4, 0, -0.4, "<-0.8"),
                         values = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(fill = "Enrichment\nScore",
         y = "Enrichment Score Rank") +
    theme(legend.position = "right",
          legend.text = element_text(size = 8),  # change legend text size
          legend.title = element_text(size = 10)) +
    theme(axis.title.x = element_blank(), # change axis title stuffs
          axis.title.y = element_text(size = 10)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1))

  if (!score_legend){  # remove legend for stacked plots
    p_strip <- p_strip + theme(legend.position = "none")
    # combine the plots for final
    p_combined <-
      cowplot::plot_grid(p_strip,
                         p_lollipop +  # remove axis and legends because it will be annotated in the strip plot
                           theme(legend.position = "none") +
                           theme(axis.text = element_blank(),
                                 axis.ticks = element_blank(),
                                 axis.title = element_blank()),

                         nrow = 1,
                         align = "h",
                         axis = "tb",
                         rel_widths = c(1, 2.5))  #note the ratio is a bit different
  } else {  # remove legend for stacked plots
    # combine the plots for final
    p_combined <-
      cowplot::plot_grid(p_lollipop +  # remove axis and legends because it will be annotated in the strip plot
                           theme(legend.position = "none") +
                           theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title = element_blank()),
                         p_strip,
                         ncol = 1,
                         align = "v",
                         axis = "lr",
                         rel_heights = c(2, 1))
  }



  return(p_combined)
}

#' Add title to combined plot
#'
#' @export
add_plot_title <- function(p_combined, plot_title, title_ratio = 10, left_margin = 0,
                           title_color = "black"){

  p_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_title,
      fontface = 'bold',
      x = 0.5,
      hjust = 0.5,
      size = 10,
      color = title_color) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, left_margin)
    )
  p_titled <- cowplot::plot_grid(p_title, p_combined, ncol = 1, rel_heights = c(1, title_ratio))
  return(p_titled)
}

