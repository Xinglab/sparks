# TODO - add DOC
calculate_mean_baseline_exp <- function(baseline_tsl){
  # calculate mean epxression data and return annotated dataframe for downstream annotation

  exp_df <- baseline_tsl@exp_df

  exp_df$gene <- unlist(lapply(rownames(exp_df), function(x) strsplit(x, "\\.")[[1]][1]))
  exp_df$info <- rownames(exp_df)

  rownames(exp_df) <- exp_df$gene
  # calculate mean exp df
  exp_df$mean <- rowMeans(exp_df[, grep("totalbaseline", colnames(exp_df))], na.rm = T)

  return(exp_df)
}

# TODO - add DOC
annotate_baseline_expression <- function(object, k562_baseline_tsl, hepg2_baseline_tsl){


  # calculate baseline mean expression
  k562_exp_df <- calculate_mean_baseline_exp(k562_baseline_tsl)
  hepg2_exp_df <- calculate_mean_baseline_exp(hepg2_baseline_tsl)

  # annotate expression
  anno_df <- object@exon_annotation$SE
  anno_df$geneID <- unlist(lapply(rownames(anno_df), function(x) strsplit(x, "\\.")[[1]][1]))

  # TODO - implement this for other SS as well

  # annotate with mean gene level expression for AS events
  anno_df$HepG2_expression <- hepg2_exp_df[anno_df$geneID, ]$mean
  anno_df$K562_expression <- k562_exp_df[anno_df$geneID, ]$mean

  object@exon_annotation$SE <- anno_df

  return(object)
}

# TODO -add DOC
subset_AS_events_with_minimal_baseline_expression_for_CLIP_analysis <- function(object, as_events, expression_threshold = 10){
  # extract available AS events and their annotation - including junction type and baseline expression
  anno_df <- object@exon_annotation$SE
  clip_events <- rownames(object@clip_result_list$SE$Upstream_5ss_intron)  # extracting any AS events with CLIP annotaetd

  # subset annotation for events tested with CLIP peak overlap
  clip_events_exp <- anno_df[clip_events,]

  # fill missing gene level expression
  clip_events_exp$HepG2_expression[is.na(clip_events_exp$HepG2_expression)] <- 0
  clip_events_exp$K562_expression[is.na(clip_events_exp$K562_expression)] <- 0

  # filter by gene level expression & known JC status
  bg_clip_events <- rownames(clip_events_exp[(clip_events_exp$HepG2_expression >= expression_threshold) &
                                               (clip_events_exp$K562_expression >= expression_threshold) &
                                               (clip_events_exp$annotation == "Known_JC"), ])

  as_events_exp_with_baseline_expression <- as_events[as_events %in% bg_clip_events]

  return(as_events_exp_with_baseline_expression)
}
