usethis::use_package("BSgenome.Hsapiens.UCSC.hg19")
usethis::use_package("VarCon")
usethis::use_package("seqinr")
usethis::use_package("ggseqlogo")

#' @export
setClass("SpliceSiteInformation",
         slots = c(ss_sequence = "list", # this would store ss seq for different event tyeps
                   ss_score = "list"))



#' Query Splice site sequence for given AS events
#'
#'
#' @param as_events AS events to get query for
#' @return splice site information object sith ss sequence and ss score
#' @export
query_splice_site_sequence_for_AS_events <- function(as_events, padding = 10){
  # load genome
  # currently genome is in GRCh37
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19


  # initialize list for data processing
  sequence_df_list = list()
  score_df_list = list()

  event_types = c("SE", "A3SS", "A5SS")
  dummy <- lapply(event_types, function(event_type){
    sequence_df_list[[event_type]] <<- list()
    score_df_list[[event_type]] <<- list()
  })

  dummy <- lapply(as_events, function(as_event){
    # debug

    # extract relevant information
    event_elements <- strsplit(as_event, ":")[[1]]
    event_type <- event_elements[11]
    strand <- event_elements[4]
    chromosome <- event_elements[3]
    if(event_type == "SE"){
      if (strand == "-"){
        upstream_5ss <- as.numeric(event_elements[9])
        cassette_3ss <- as.numeric(event_elements[6])
        cassette_5ss <- as.numeric(event_elements[5])
        dnstream_3ss <- as.numeric(event_elements[8])

        # retreive sequence
        upst_5ss_seq <- reverseComplement(genome[[chromosome]][(upstream_5ss - padding + 1) : (upstream_5ss + padding)])
        cass_3ss_seq <- reverseComplement(genome[[chromosome]][(cassette_3ss - padding + 1) : (cassette_3ss + padding)])
        cass_5ss_seq <- reverseComplement(genome[[chromosome]][(cassette_5ss - padding + 1) : (cassette_5ss + padding)])
        dnst_3ss_seq <- reverseComplement(genome[[chromosome]][(dnstream_3ss - padding + 1) : (dnstream_3ss + padding)])

        # calculate score - 3E/6I for 3'ss, 20I/3E for 5'ss
        upst_5ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(upstream_5ss - 5) : (upstream_5ss + 3)]), 5)
        cass_3ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(cassette_3ss - 2) : (cassette_3ss + 20)]), 3)
        cass_5ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(cassette_5ss - 5) : (cassette_5ss + 3)]), 5)
        dnst_3ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(dnstream_3ss - 2) : (dnstream_3ss + 20)]), 3)



      } else if (strand == "+"){
        upstream_5ss <- as.numeric(event_elements[8])
        cassette_3ss <- as.numeric(event_elements[5])
        cassette_5ss <- as.numeric(event_elements[6])
        dnstream_3ss <- as.numeric(event_elements[9])

        upst_5ss_seq <- genome[[chromosome]][(upstream_5ss - padding + 1) : (upstream_5ss + padding)]
        cass_3ss_seq <- genome[[chromosome]][(cassette_3ss - padding + 1) : (cassette_3ss + padding)]
        cass_5ss_seq <- genome[[chromosome]][(cassette_5ss - padding + 1) : (cassette_5ss + padding)]
        dnst_3ss_seq <- genome[[chromosome]][(dnstream_3ss - padding + 1) : (dnstream_3ss + padding)]

        # calculate score - 3E/6I for 3'ss, 20I/3E for 5'ss
        upst_5ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(upstream_5ss - 2) : (upstream_5ss + 6)], 5)
        cass_3ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(cassette_3ss - 19) : (cassette_3ss + 3)], 3)
        cass_5ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(cassette_5ss - 2) : (cassette_5ss + 6)], 5)
        dnst_3ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(dnstream_3ss - 19) : (dnstream_3ss + 3)], 3)
      }
      # aggergate sequence information
      sequence_df <- data.frame(upstream_5ss_seq = as.character(upst_5ss_seq),
                                cassette_3ss_seq = as.character(cass_3ss_seq),
                                cassette_5ss_seq = as.character(cass_5ss_seq),
                                dnstream_3ss_seq = as.character(dnst_3ss_seq))
      score_df <- data.frame(upstream_5ss_seq = as.numeric(upst_5ss_maxent_score),
                             cassette_3ss_seq = as.numeric(cass_3ss_maxent_score),
                             cassette_5ss_seq = as.numeric(cass_5ss_maxent_score),
                             dnstream_3ss_seq = as.numeric(dnst_3ss_maxent_score))
      rownames(sequence_df) <- as_event
      rownames(score_df) <- as_event

      # add them to outside list for calculation
      sequence_df_list[[event_type]][[length(sequence_df_list[[event_type]]) + 1]] <<- sequence_df
      score_df_list[[event_type]][[length(score_df_list[[event_type]]) + 1]] <<- score_df

    } else if (event_type == "A3SS"){
      if (strand == "+"){
        # example = "ENSG00000126456.15_3:IRF3:chr19:-:50166461:50166771:50166461:50166515:50168887:50168962:A3SS"

        # define splice sites
        upstream_5ss <- as.numeric(event_elements[10])
        long_3ss <- as.numeric(event_elements[5])
        short_3ss <- as.numeric(event_elements[7])

        # retreive sequence
        upst_5ss_seq <- genome[[chromosome]][(upstream_5ss - padding + 1) : (upstream_5ss + padding)]
        long_3ss_seq <- genome[[chromosome]][(long_3ss - padding + 1) : (long_3ss + padding)]
        shrt_3ss_seq <- genome[[chromosome]][(short_3ss - padding + 1) : (short_3ss + padding)]

        # calculate score - 3E/6I for 3'ss, 20I/3E for 5'ss
        upst_5ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(upstream_5ss - 2) : (upstream_5ss + 6)], 5)
        long_3ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(long_3ss - 19) : (long_3ss + 3)], 3)
        shrt_3ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(short_3ss - 19) : (short_3ss + 3)], 3)


      } else if (strand == "-"){
        # define splice sites
        upstream_5ss <- as.numeric(event_elements[9])
        long_3ss <- as.numeric(event_elements[6])
        short_3ss <- as.numeric(event_elements[8])

        # retreive sequence
        upst_5ss_seq <- reverseComplement(genome[[chromosome]][(upstream_5ss - padding + 1) : (upstream_5ss + padding)])
        long_3ss_seq <- reverseComplement(genome[[chromosome]][(long_3ss - padding + 1) : (long_3ss + padding)])
        shrt_3ss_seq <- reverseComplement(genome[[chromosome]][(short_3ss - padding + 1) : (short_3ss + padding)])

        # calculate score - 3E/6I for 3'ss, 20I/3E for 5'ss
        upst_5ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(upstream_5ss - 5) : (upstream_5ss + 3)]), 5)
        long_3ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(long_3ss - 2) : (long_3ss + 20)]), 3)
        shrt_3ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(short_3ss - 2) : (short_3ss + 20)]), 3)
      }
      # aggergate sequence information
      sequence_df <- data.frame(upstream_5ss_seq = as.character(upst_5ss_seq),
                                long_3ss_seq = as.character(long_3ss_seq),
                                short_3ss_seq = as.character(shrt_3ss_seq))
      score_df <- data.frame(upstream_5ss_seq = as.numeric(upst_5ss_maxent_score),
                             long_3ss_seq = as.numeric(long_3ss_maxent_score),
                             short_3ss_seq = as.numeric(shrt_3ss_maxent_score))
      rownames(sequence_df) <- as_event
      rownames(score_df) <- as_event
      sequence_df_list[[event_type]][[length(sequence_df_list[[event_type]]) + 1]] <<- sequence_df
      score_df_list[[event_type]][[length(score_df_list[[event_type]]) + 1]] <<- score_df

    }
    else if (event_type == "A5SS") {
      if (strand == "+"){
        # example = "ENSG00000108107.14_3:RPL28:chr19:+:55899364:55899852:55899364:55899416:55902921:55903171:A5SS"

        # define splice sites
        long_5ss <- as.numeric(event_elements[6])
        short_5ss <- as.numeric(event_elements[8])
        dnstream_3ss <- as.numeric(event_elements[9])

        # retreive sequence
        long_5ss_seq <- genome[[chromosome]][(long_5ss - padding + 1) : (long_5ss + padding)]
        shrt_5ss_seq <- genome[[chromosome]][(short_5ss - padding + 1) : (short_5ss + padding)]
        dnst_3ss_seq <- genome[[chromosome]][(dnstream_3ss - padding + 1) : (dnstream_3ss + padding)]

        # calculate score - 3E/6I for 3'ss, 20I/3E for 5'ss
        long_5ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(long_5ss - 2) : (long_5ss + 6)], 5)
        shrt_5ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(short_5ss - 2) : (short_5ss + 6)], 5)
        dnst_3ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(dnstream_3ss - 19) : (dnstream_3ss + 3)], 3)


      } else if (strand == "-"){
        # example = "ENSG00000142089.15_3:IFITM3:chr11:-:320464:320683:320564:320683:320112:320248:A5SS"
        # define splice sites
        long_5ss <- as.numeric(event_elements[5])
        short_5ss <- as.numeric(event_elements[7])
        dnstream_3ss <- as.numeric(event_elements[10])

        # retreive sequence
        long_3ss_seq <- reverseComplement(genome[[chromosome]][(long_5ss - padding + 1) : (long_5ss + padding)])
        shrt_3ss_seq <- reverseComplement(genome[[chromosome]][(short_5ss - padding + 1) : (short_5ss + padding)])
        dnst_5ss_seq <- reverseComplement(genome[[chromosome]][(dnstream_3ss - padding + 1) : (dnstream_3ss + padding)])

        # calculate score - 3E/6I for 3'ss, 20I/3E for 5'ss
        long_5ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(long_5ss - 5) : (long_5ss + 3)]), 5)
        shrt_5ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(short_5ss - 5) : (short_5ss + 3)]), 5)
        dnst_3ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(dnstream_3ss - 2) : (dnstream_3ss + 20)]), 3)
      }
      # aggergate sequence information
      sequence_df <- data.frame(long_5ss_seq = as.character(long_5ss_seq),
                                short_5ss_seq = as.character(shrt_5ss_seq),
                                dnstream_3ss_seq = as.character(dnst_3ss_seq))
      score_df <- data.frame(long_5ss_seq = as.numeric(long_5ss_maxent_score),
                             short_5ss_seq = as.numeric(shrt_5ss_maxent_score),
                             dnstream_3ss_seq = as.numeric(dnst_3ss_maxent_score))
      rownames(sequence_df) <- as_event
      rownames(score_df) <- as_event
      sequence_df_list[[event_type]][[length(sequence_df_list[[event_type]]) + 1]] <<- sequence_df
      score_df_list[[event_type]][[length(score_df_list[[event_type]]) + 1]] <<- score_df
    }})
  # store them in new splice site information object
  object <- new("SpliceSiteInformation")
  dummy <- lapply(event_types, function(event_type){
    if (length(sequence_df_list[[event_type]]) > 0){ # if the list is not empty - possible for specific events
      # assuming score df is the same - as it should be
      print(event_type)

      # combine data for each event type
      ss_sequence_df <- do.call(rbind, sequence_df_list[[event_type]])
      ss_score_df <- do.call(rbind, score_df_list[[event_type]])
      # add them to the object
      object@ss_sequence[[event_type]] <<- ss_sequence_df
      object@ss_score[[event_type]] <<- ss_score_df

    }
  })
  return(object)
}



#' Generate Splice site sequence logo
#'
#'
#' @param object splice site information object with ss sequence and ss score
#' @return list of plots with sequence logo
#' @export
generate_ss_sequence_logo <- function(object, event_type, padding = 10){
  # calcalate mean scores
  ss_scores <- colMeans(object@ss_score[[event_type]])
  # break here? TODO
  # generate consensus motif plot
  ss_sequence_df <- object@ss_sequence[[event_type]]
  ss_events <- colnames(ss_sequence_df)
  p_list <- lapply(ss_events, function(sequence_name){

    dummy <- ss_sequence_df[, sequence_name, drop=FALSE]
    mat <- as.data.frame(t(do.call(cbind, apply(dummy, 2, function(x) as.data.frame(strsplit(x, split = ""))))))

    consensus_pwm <- consensus(as.matrix(mat), method = "profile")

    # determine splice site type
    # TODO - fix for other A3SS/A5SS
    # if (strsplit(sequence_name, "_")[[1]][2] == '5ss'){
    #   seq_label <- c(seq(10, 1), seq(-1, -10))
    # } else if (strsplit(sequence_name, "_")[[1]][2] == '3ss'){
    #   seq_label <- c(seq(-10, -1), seq(1, 10))
    # }
    seq_label <- c(seq(-padding, -1), seq(1, padding))
    # determine mean score
    mean_score_text <- sprintf("Mean MaxEntScan Score = %.*f", 3, round(ss_scores[sequence_name], digits = 3))

    p <- ggseqlogo(consensus_pwm) +
      scale_x_continuous(labels = seq_label, breaks = seq(padding * 2)) +
      geom_vline(xintercept = padding + 0.5, color = "red", linetype = "dashed") +
      labs(subtitle = mean_score_text)+
      theme(plot.subtitle = element_text(size = 12))
    return(p)
  })

  return(p_list)
}



#' Generate Splice Site Strength Density plot
#'
#'
#' @param inc_object splice site information object with ss sequence and ss score for inclusion part (TOP)
#' @param exc_object splice site information object with ss sequence and ss score for inclusion part (BOTTOM)
#' @param event_type Event Type - SE, A3SS, A5SS
#' @return list of plots with sequence strength distribution
#' @export
generate_ss_sequence_score <- function(inc_object, exc_object, event_type){
  # calcalate mean scores
  inc_ss_score_df <- inc_object@ss_score[[event_type]]
  exc_ss_score_df <- exc_object@ss_score[[event_type]]
  # generate consensus motif plot
  ss_sequence_df <- inc_object@ss_sequence[[event_type]]
  ss_events <- colnames(ss_sequence_df)
  p_list <- lapply(ss_events, function(sequence_name){

    inc_score <- inc_ss_score_df[, sequence_name, drop=FALSE]
    exc_score <- exc_ss_score_df[, sequence_name, drop=FALSE]

    # add anntation
    inc_score$direction <- "inclusion"
    exc_score$direction <- "exclusion"

    score_df <- rbind(inc_score, exc_score)
    colnames(score_df) <- c("score", "psi")

    p <- ggplot(score_df, aes(x = score, y = psi, fill = psi)) +
      geom_violin() +
      geom_boxplot()

    return(p)
  })

  return(p_list)
}



#' Combining script to generate Splice Site Analysis Plot
#'
#'
#' @param p_list_inc list of plots for positive delta psi
#' @param p_list_exc list of plots for negative delta psi
#' @param p_list_score list of plots for splice site score distribution
#' @param condition_label labels to be added in the plot
#' @return combined plot with all the plots merged
#' @export
generate_combined_splice_site_plot <- function(p_list_inc, p_list_exc, p_list_score, condition_label_one, condition_label_two){
  p <- cowplot::plot_grid(p_list_inc[[1]] + labs(y = sprintf("%s\nBits", condition_label_one)),
                 p_list_inc[[2]] + theme(axis.title.y = element_blank()),
                 p_list_inc[[3]] + theme(axis.title.y = element_blank()),
                 p_list_inc[[4]] + theme(axis.title.y = element_blank()),
                 p_list_score[[1]] +
                   labs(y = "MaxEntScan Score") +
                   theme(legend.position = "none"),
                 p_list_score[[2]] + theme(axis.title.y = element_blank()) +
                   theme(legend.position = "none") +
                   theme(axis.text.y = element_blank()),
                 p_list_score[[3]] + theme(axis.title.y = element_blank()) +
                   theme(legend.position = "none")+
                   theme(axis.text.y = element_blank()),
                 p_list_score[[4]] + theme(axis.title.y = element_blank()) +
                   theme(legend.position = "none")+
                   theme(axis.text.y = element_blank()),
                 p_list_exc[[1]] +
                   labs(y = sprintf("%s\nBits", condition_label_two)) +
                   labs(x = "Upstream 5'SS"),
                 p_list_exc[[2]] + theme(axis.title.y = element_blank()) +
                   labs(x = "Cassette 3'SS"),
                 p_list_exc[[3]] + theme(axis.title.y = element_blank()) +
                   labs(x = "Cassette 5'SS"),
                 p_list_exc[[4]] + theme(axis.title.y = element_blank()) +
                   labs(x = "Downstream 3'SS"),

                 nrow = 3)
  return(p)
}



#' Main Wrapper to generate Splice Site Analysis Plot
#'
#'
#' @param pos_events AS events with positive delta PSI - this will be on the top when plots are generated
#' @param neg_events AS events with negative delta PSI - this will be on the bot when plots are generated
#' @param condition_label labels to be added in the plot
#' @param event_type Event Type - SE, A3SS, A5SS
#' @return combined plot with all the plots merged
#' @export
generate_splice_site_consensus_plot <- function(pos_events, neg_events, condition_label_one, condition_label_two, event_type = "SE", padding = 10){

  inc_ss <- query_splice_site_sequence_for_AS_events(pos_events, padding = padding)
  exc_ss <- query_splice_site_sequence_for_AS_events(neg_events, padding = padding)
  p_list_inc <- generate_ss_sequence_logo(inc_ss, event_type, padding = padding)
  p_list_exc <- generate_ss_sequence_logo(exc_ss, event_type, padding = padding)
  p_list_score <- generate_ss_sequence_score(inc_ss, exc_ss, event_type)

  combined_p <- generate_combined_splice_site_plot(p_list_inc, p_list_exc, p_list_score, condition_label_one, condition_label_two)
  return(combined_p)
}




##### DEPRECATED FUNCTIONS ####
query_splice_site_sequence_for_AS_events_four_field <- function(as_events){

  padding <- 10

  # initialize list for data processing
  sequence_df_list = list()
  score_df_list = list()

  event_types = c("SE")
  dummy <- lapply(event_types, function(event_type){
    sequence_df_list[[event_type]] <<- list()
    score_df_list[[event_type]] <<- list()
  })

  dummy <- lapply(as_events, function(as_event){
    # debug

    # extract relevant information
    event_elements <- strsplit(as_event, "_")[[1]]
    # event_type <- event_elements[11]
    event_type <- "SE" # seems like PNAS paper only includes SE events and it is not annotated
    strand <- event_elements[3]
    chromosome <- event_elements[2]

    if(event_type == "SE"){
      if (strand == "-"){
        upstream_5ss <- as.numeric(event_elements[7])
        cassette_3ss <- as.numeric(event_elements[5])
        cassette_5ss <- as.numeric(event_elements[4])
        dnstream_3ss <- as.numeric(event_elements[6])

        # retreive sequence
        upst_5ss_seq <- reverseComplement(genome[[chromosome]][(upstream_5ss - padding + 1) : (upstream_5ss + padding)])
        cass_3ss_seq <- reverseComplement(genome[[chromosome]][(cassette_3ss - padding + 1) : (cassette_3ss + padding)])
        cass_5ss_seq <- reverseComplement(genome[[chromosome]][(cassette_5ss - padding + 1) : (cassette_5ss + padding)])
        dnst_3ss_seq <- reverseComplement(genome[[chromosome]][(dnstream_3ss - padding + 1) : (dnstream_3ss + padding)])

        # calculate score - 3E/6I for 3'ss, 20I/3E for 5'ss
        upst_5ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(upstream_5ss - 5) : (upstream_5ss + 3)]), 5)
        cass_3ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(cassette_3ss - 2) : (cassette_3ss + 20)]), 3)
        cass_5ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(cassette_5ss - 5) : (cassette_5ss + 3)]), 5)
        dnst_3ss_maxent_score <- calculateMaxEntScanScore(reverseComplement(genome[[chromosome]][(dnstream_3ss - 2) : (dnstream_3ss + 20)]), 3)



      } else if (strand == "+"){
        upstream_5ss <- as.numeric(event_elements[6])
        cassette_3ss <- as.numeric(event_elements[4])
        cassette_5ss <- as.numeric(event_elements[5])
        dnstream_3ss <- as.numeric(event_elements[7])

        upst_5ss_seq <- genome[[chromosome]][(upstream_5ss - padding + 1) : (upstream_5ss + padding)]
        cass_3ss_seq <- genome[[chromosome]][(cassette_3ss - padding + 1) : (cassette_3ss + padding)]
        cass_5ss_seq <- genome[[chromosome]][(cassette_5ss - padding + 1) : (cassette_5ss + padding)]
        dnst_3ss_seq <- genome[[chromosome]][(dnstream_3ss - padding + 1) : (dnstream_3ss + padding)]

        # calculate score - 3E/6I for 3'ss, 20I/3E for 5'ss
        upst_5ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(upstream_5ss - 2) : (upstream_5ss + 6)], 5)
        cass_3ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(cassette_3ss - 19) : (cassette_3ss + 3)], 3)
        cass_5ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(cassette_5ss - 2) : (cassette_5ss + 6)], 5)
        dnst_3ss_maxent_score <- calculateMaxEntScanScore(genome[[chromosome]][(dnstream_3ss - 19) : (dnstream_3ss + 3)], 3)
      }
      # aggergate sequence information
      sequence_df <- data.frame(upstream_5ss_seq = as.character(upst_5ss_seq),
                                cassette_3ss_seq = as.character(cass_3ss_seq),
                                cassette_5ss_seq = as.character(cass_5ss_seq),
                                dnstream_3ss_seq = as.character(dnst_3ss_seq))
      score_df <- data.frame(upstream_5ss_seq = as.numeric(upst_5ss_maxent_score),
                             cassette_3ss_seq = as.numeric(cass_3ss_maxent_score),
                             cassette_5ss_seq = as.numeric(cass_5ss_maxent_score),
                             dnstream_3ss_seq = as.numeric(dnst_3ss_maxent_score))
      rownames(sequence_df) <- as_event
      rownames(score_df) <- as_event

      # add them to outside list for calculation
      sequence_df_list[[event_type]][[length(sequence_df_list[[event_type]]) + 1]] <<- sequence_df
      score_df_list[[event_type]][[length(score_df_list[[event_type]]) + 1]] <<- score_df

    }})
  # store them in new splice site information object
  object <- new("SpliceSiteInformation")
  dummy <- lapply(event_types, function(event_type){
    if (length(sequence_df_list[[event_type]]) > 0){ # if the list is not empty - possible for specific events
      # assuming score df is the same - as it should be
      print(event_type)

      # combine data for each event type
      ss_sequence_df <- do.call(rbind, sequence_df_list[[event_type]])
      ss_score_df <- do.call(rbind, score_df_list[[event_type]])
      # add them to the object
      object@ss_sequence[[event_type]] <<- ss_sequence_df
      object@ss_score[[event_type]] <<- ss_score_df

    }
  })
  return(object)
}



generate_ss_sequence_logo_four_field <- function(object, event_type){
  # calcalate mean scores
  ss_scores <- colMeans(object@ss_score[[event_type]])
  # break here? TODO
  # generate consensus motif plot
  ss_sequence_df <- object@ss_sequence[[event_type]]
  ss_events <- colnames(ss_sequence_df)
  p_list <- lapply(ss_events, function(sequence_name){

    dummy <- ss_sequence_df[, sequence_name, drop=FALSE]
    mat <- as.data.frame(t(do.call(cbind, apply(dummy, 2, function(x) as.data.frame(strsplit(x, split = ""))))))

    consensus_pwm <- consensus(as.matrix(mat), method = "profile")

    # determine splice site type
    # TODO - fix for other A3SS/A5SS
    if (strsplit(sequence_name, "_")[[1]][2] == '5ss'){
      seq_label <- c(seq(10, 1), seq(-1, -10))
    } else if (strsplit(sequence_name, "_")[[1]][2] == '3ss'){
      seq_label <- c(seq(-10, -1), seq(1, 10))
    }

    # determine mean score
    mean_score_text <- sprintf("Mean MaxEntScan Score = %f", ss_scores[sequence_name])

    p <- ggseqlogo(consensus_pwm) +
      scale_x_continuous(labels = seq_label, breaks = seq(20)) +
      geom_vline(xintercept = 10.5, color = "red", linetype = "dashed") +
      labs(subtitle = mean_score_text)+
      theme(plot.subtitle = element_text(size = 8))
    return(p)
  })

  return(p_list)
}




# TODO - add doc
calculate_ss_score_diff_directional <- function(pos_events, neg_events){

  inc_ss <- query_splice_site_sequence_for_AS_events(pos_events)
  exc_ss <- query_splice_site_sequence_for_AS_events(neg_events)
  # calculate differential scores
  score_diff_5ss_inc <- as.data.frame(inc_ss@ss_score$SE$cassette_5ss_seq - inc_ss@ss_score$SE$upstream_5ss_seq)
  score_diff_3ss_inc <- as.data.frame(inc_ss@ss_score$SE$cassette_3ss_seq - inc_ss@ss_score$SE$dnstream_3ss_seq)
  score_diff_5ss_exc <- as.data.frame(exc_ss@ss_score$SE$cassette_5ss_seq - exc_ss@ss_score$SE$upstream_5ss_seq)
  score_diff_3ss_exc <- as.data.frame(exc_ss@ss_score$SE$cassette_3ss_seq - exc_ss@ss_score$SE$dnstream_3ss_seq)

  # fix the column name
  colnames(score_diff_5ss_inc) <- "differential"
  colnames(score_diff_3ss_inc) <- "differential"
  colnames(score_diff_5ss_exc) <- "differential"
  colnames(score_diff_3ss_exc) <- "differential"

  # calculate average
  score_diff_5ss_inc$average <- (inc_ss@ss_score$SE$upstream_5ss_seq + inc_ss@ss_score$SE$cassette_5ss_seq)/2
  score_diff_3ss_inc$average <- (inc_ss@ss_score$SE$cassette_3ss_seq + inc_ss@ss_score$SE$dnstream_3ss_seq)/2
  score_diff_5ss_exc$average <- (exc_ss@ss_score$SE$upstream_5ss_seq + exc_ss@ss_score$SE$cassette_5ss_seq)/2
  score_diff_3ss_exc$average <- (exc_ss@ss_score$SE$cassette_3ss_seq + exc_ss@ss_score$SE$dnstream_3ss_seq)/2

  # annotate direction
  score_diff_5ss_inc$direction <- "positive"
  score_diff_3ss_inc$direction <- "positive"
  score_diff_5ss_exc$direction <- "negative"
  score_diff_3ss_exc$direction <- "negative"

  # annotate splice site
  score_diff_5ss_inc$site <- '5ss'
  score_diff_3ss_inc$site <- '3ss'
  score_diff_5ss_exc$site <- '5ss'
  score_diff_3ss_exc$site <- '3ss'

  # annotate events
  score_diff_5ss_inc$event <- pos_events
  score_diff_3ss_inc$event <- pos_events
  score_diff_5ss_exc$event <- neg_events
  score_diff_3ss_exc$event <- neg_events

  # combine everything
  score_diff_plot_df <- rbind(score_diff_5ss_inc,
                              score_diff_3ss_inc,
                              score_diff_5ss_exc,
                              score_diff_3ss_exc)

  return(score_diff_plot_df)
}




#' Query Splice site sequence for given AS events
#'
#'
#' @param as_events AS events to get query for
#' @return splice site information object sith ss sequence and ss score
#' @export
calculate_polypyrimidine_base_frequency <- function(as_events, padding = 100){
  # load genome
  # currently genome is in GRCh37
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19


  # initialize list for data processing
  event_types = c("SE", "A3SS", "A5SS")

  pyr_result_list <- lapply(as_events, function(as_event){
    # debug

    # extract relevant information
    event_elements <- strsplit(as_event, ":")[[1]]
    event_type <- event_elements[11]
    strand <- event_elements[4]
    chromosome <- event_elements[3]
    if(event_type == "SE"){
      if (strand == "-"){
        upstream_5ss <- as.numeric(event_elements[9])
        cassette_3ss <- as.numeric(event_elements[6])
        cassette_5ss <- as.numeric(event_elements[5])
        dnstream_3ss <- as.numeric(event_elements[8])

        # retreive sequence
        cass_3ss_seq <- reverseComplement(genome[[chromosome]][(cassette_3ss) : (cassette_3ss + padding)])
        dnst_3ss_seq <- reverseComplement(genome[[chromosome]][(dnstream_3ss) : (dnstream_3ss + padding)])

        # calculate frequency
        cass_3ss_freq <- alphabetFrequency(cass_3ss_seq)[c("A","C","G","T")]/sum(alphabetFrequency(cass_3ss_seq)[c("A","C","G","T")])
        dnst_3ss_freq <- alphabetFrequency(dnst_3ss_seq)[c("A","C","G","T")]/sum(alphabetFrequency(dnst_3ss_seq)[c("A","C","G","T")])

      } else if (strand == "+"){
        upstream_5ss <- as.numeric(event_elements[8])
        cassette_3ss <- as.numeric(event_elements[5])
        cassette_5ss <- as.numeric(event_elements[6])
        dnstream_3ss <- as.numeric(event_elements[9])

        cass_3ss_seq <- genome[[chromosome]][(cassette_3ss - padding) : (cassette_3ss)]
        dnst_3ss_seq <- genome[[chromosome]][(dnstream_3ss - padding) : (dnstream_3ss)]

        # calculate frequency
        cass_3ss_freq <- alphabetFrequency(cass_3ss_seq)[c("A","C","G","T")]/sum(alphabetFrequency(cass_3ss_seq)[c("A","C","G","T")])
        dnst_3ss_freq <- alphabetFrequency(dnst_3ss_seq)[c("A","C","G","T")]/sum(alphabetFrequency(dnst_3ss_seq)[c("A","C","G","T")])

      }
      ## aggregate frequecy information
      # put it in dataframe for processing
      cass_3ss_df <- as.data.frame(cass_3ss_freq)
      dnst_3ss_df <- data.frame(dnst_3ss_freq)
      colnames(cass_3ss_df) <- c("freq")
      colnames(dnst_3ss_df) <- c("freq")

      # annotate
      cass_3ss_df$base <- rownames(cass_3ss_df)
      dnst_3ss_df$base <- rownames(dnst_3ss_df)
      rownames(cass_3ss_df) <- NULL
      rownames(dnst_3ss_df) <- NULL

      cass_3ss_df$region <- "cassette_3ss_Pyr"
      dnst_3ss_df$region <- "dnstream_3ss_Pyr"

      # bind them together
      pyr_freq_df <- rbind(cass_3ss_df, dnst_3ss_df)
      pyr_freq_df$event <- as_event
      return(pyr_freq_df)
    }})
  pyr_result_df <- do.call(rbind, pyr_result_list)
  return(pyr_result_df)
}



#' Query Splice site sequence for given AS events
#'
#'
#' @param as_events AS events to get query for
#' @return splice site information object sith ss sequence and ss score
#' @export
calculate_motif_occurrence <- function(as_events, padding = 100, motif = "CCCCC"){

  # load genome
  # currently genome is in GRCh37
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19


  # initialize list for data processing
  event_types = c("SE", "A3SS", "A5SS")

  pyr_result_list <- lapply(as_events, function(as_event){
    # debug

    # extract relevant information
    event_elements <- strsplit(as_event, ":")[[1]]
    event_type <- event_elements[11]
    strand <- event_elements[4]
    chromosome <- event_elements[3]
    if(event_type == "SE"){
      if (strand == "-"){
        upstream_5ss <- as.numeric(event_elements[9])
        cassette_3ss <- as.numeric(event_elements[6])
        cassette_5ss <- as.numeric(event_elements[5])
        dnstream_3ss <- as.numeric(event_elements[8])

        # retreive sequence
        cass_3ss_seq <- reverseComplement(genome[[chromosome]][(cassette_3ss) : (cassette_3ss + padding)])
        dnst_3ss_seq <- reverseComplement(genome[[chromosome]][(dnstream_3ss) : (dnstream_3ss + padding)])

      } else if (strand == "+"){
        upstream_5ss <- as.numeric(event_elements[8])
        cassette_3ss <- as.numeric(event_elements[5])
        cassette_5ss <- as.numeric(event_elements[6])
        dnstream_3ss <- as.numeric(event_elements[9])

        cass_3ss_seq <- genome[[chromosome]][(cassette_3ss - padding) : (cassette_3ss)]
        dnst_3ss_seq <- genome[[chromosome]][(dnstream_3ss - padding) : (dnstream_3ss)]

      }
      # aggregate frequecy information
      cass_3ss_count <- countPattern(motif, cass_3ss_seq)
      dnst_3ss_count <- countPattern(motif, dnst_3ss_seq)

      result_df <- data.frame(count = c(cass_3ss_count, dnst_3ss_count),
                              region = c("cassette_3ss_Pyr", "dnstream_3ss_Pyr"))
      result_df$event <- as_event
      result_df$motif <- motif
      return(result_df)
    }})
  pyr_result_df <- do.call(rbind, pyr_result_list)
  return(pyr_result_df)
}
