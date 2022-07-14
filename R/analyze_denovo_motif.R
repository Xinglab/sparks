
##### STREME method ######

run_directional_STREME <- function(pos_events, neg_events, output_dir){
  # make output dir if doesnt exist
  if (!dir.exists(output_dir)) dir.create(output_dir)

  # run enrichment test for all regions seperately
  region_list <- c("upstream_5ss_exon",
                   "upstream_5ss_intron",
                   "cassette_3ss_intron",
                   "cassette_3ss_exon",
                   "cassette_5ss_exon",
                   "cassette_5ss_intron",
                   "dnstream_3ss_intron",
                   "dnstream_3ss_exon")
  region_result_list <- lapply(region_list, function(region){
    streme_output_dir <- sprintf("%s/%s/", output_dir, region)
    if (!dir.exists(streme_output_dir)) dir.create(streme_output_dir)

    # query fasta sequence
    print("Fetching Sequences")
    pos_event_fasta <- do.call(rbind, lapply(pos_events, function(x) call_fasta_from_event(x, region)))
    neg_event_fasta <- do.call(rbind, lapply(neg_events, function(x) call_fasta_from_event(x, region)))

    # write csv files
    print("Writing Sequences")
    pos_seq_fasta <- sprintf("%s/positive_events.fasta", streme_output_dir)
    neg_seq_fasta <- sprintf("%s/negative_events.fasta", streme_output_dir)

    write.table(pos_event_fasta, pos_seq_fasta, sep = '\n', quote = F, row.names = F, col.names = F)
    write.table(neg_event_fasta, neg_seq_fasta, sep = '\n', quote = F, row.names = F, col.names = F)


    # run streme
    print("Running STREME")
    pos_output_dir <- sprintf("%s/positive/", streme_output_dir)
    neg_output_dir <- sprintf("%s/negative/", streme_output_dir)


    cmd <- sprintf("~/meme/bin/streme --p %s --n %s --objfun de --oc %s --rna --minw 5 --maxw 10 --patience 10", pos_seq_fasta, neg_seq_fasta, pos_output_dir)
    system(cmd)
    cmd <- sprintf("~/meme/bin/streme --p %s --n %s --objfun de --oc %s --rna --minw 5 --maxw 10 --patience 10", neg_seq_fasta, pos_seq_fasta, neg_output_dir)
    system(cmd)

    # read in the result
    print("Gathering STREME results")
    pos_result_xml_file <- sprintf("%s/streme.xml", pos_output_dir)
    neg_result_xml_file <- sprintf("%s/streme.xml", neg_output_dir)

    pos_streme_pval_df <- import_streme_motif_pval_df(pos_result_xml_file)
    neg_streme_pval_df <- import_streme_motif_pval_df(neg_result_xml_file)

    # annotate
    pos_streme_pval_df$direction <- "positive"
    pos_streme_pval_df$region <- region
    neg_streme_pval_df$direction <- "negative"
    neg_streme_pval_df$region <- region


    if (is.null(dim(pos_streme_pval_df)) & is.null(dim(neg_streme_pval_df))){
      return(NULL)
    } else if (is.null(dim(pos_streme_pval_df))) {
      return(neg_streme_pval_df)
    } else if (is.null(dim(neg_streme_pval_df))) {
      return(pos_streme_pval_df)
    } else {
      streme_pval_df <- rbind(pos_streme_pval_df, neg_streme_pval_df)
      return(streme_pval_df)
    }

  })

  region_result_df <- do.call(rbind, region_result_list)
  return(region_result_df)
}

import_streme_motif_pval_df <- function(xml_result_file, pval_threshold = 0.05){
  xml_result <- xmlToList(xmlParse(xml_result_file))
  sig_motif_list <- list()

  dummy <- lapply(xml_result$motifs, function(motif_entry){


    # extract pwm
    motif_len <- as.numeric(motif_entry$.attrs["width"])

    motif_pwm <- as.data.frame(do.call(rbind, lapply(seq(motif_len), function(position){
      motif_pwm_pos <- as.numeric(motif_entry[[position]])
      return(motif_pwm_pos)
    })))
    colnames(motif_pwm) <- c("A", "C", "G", "U")

    motif_pwm_logo <- ggseqlogo(t(motif_pwm)) # TODO - currently this is not utilized anywhere

    # extract pval and motif
    motif_sequence <- strsplit(motif_entry$.attrs["id"], "-")[[1]][2]
    motif_pval <- as.numeric(motif_entry$.attrs["test_pvalue"])

    if (motif_pval < pval_threshold){
      motif_df <- data.frame(sequence = motif_sequence,
                             pval = motif_pval)
      sig_motif_list[[length(sig_motif_list) + 1 ]] <<- motif_df
      return()
    }
  })

  motif_pval_df <-  do.call(rbind, sig_motif_list)
  return(motif_pval_df)
}

# extract fastq sequence
call_fasta_from_event <- function(event, event_type = "exon",
                                  exon_padding = 50,
                                  intron_padding = 300){

  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  bed_region_list <- strsplit(event, ":")[[1]]
  chr <- bed_region_list[3]
  direction <- bed_region_list[4]
  if (event_type == "upstream_5ss_exon"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[9]

      fasta_seq <- reverseComplement(genome[[chr]][(as.numeric(splice_site_coord)) : (as.numeric(splice_site_coord) + exon_padding)])
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[8]

      fasta_seq <- genome[[chr]][(as.numeric(splice_site_coord) - exon_padding) : (as.numeric(splice_site_coord))]
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))

    }
  } else if (event_type == "upstream_5ss_intron"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[9]

      fasta_seq <- reverseComplement(genome[[chr]][(as.numeric(splice_site_coord) - intron_padding) : (as.numeric(splice_site_coord))])
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))
    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[8]

      fasta_seq <- genome[[chr]][(as.numeric(splice_site_coord)) : (as.numeric(splice_site_coord) + intron_padding)]
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))
    }
  } else if (event_type == "cassette_3ss_intron"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[6]

      fasta_seq <- reverseComplement(genome[[chr]][(as.numeric(splice_site_coord)) : (as.numeric(splice_site_coord) + intron_padding)])
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[5]

      fasta_seq <- genome[[chr]][(as.numeric(splice_site_coord) - intron_padding) : (as.numeric(splice_site_coord))]
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))

    }
  } else if (event_type == "cassette_3ss_exon"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[6]

      fasta_seq <- reverseComplement(genome[[chr]][(as.numeric(splice_site_coord) - exon_padding) : (as.numeric(splice_site_coord))])
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[5]

      fasta_seq <- genome[[chr]][(as.numeric(splice_site_coord)) : (as.numeric(splice_site_coord) + exon_padding)]
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))

    }
  } else if (event_type == "cassette_5ss_exon"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[5]

      fasta_seq <- reverseComplement(genome[[chr]][(as.numeric(splice_site_coord)) : (as.numeric(splice_site_coord) + exon_padding)])
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[6]

      fasta_seq <- genome[[chr]][(as.numeric(splice_site_coord) - exon_padding) : (as.numeric(splice_site_coord))]
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))
    }
  } else if (event_type == "cassette_5ss_intron"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[5]

      fasta_seq <- reverseComplement(genome[[chr]][(as.numeric(splice_site_coord) - intron_padding) : (as.numeric(splice_site_coord))])
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))
    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[6]

      fasta_seq <- genome[[chr]][(as.numeric(splice_site_coord)) : (as.numeric(splice_site_coord) + intron_padding)]
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))
    }
  } else if (event_type == "dnstream_3ss_intron"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[8]

      fasta_seq <- reverseComplement(genome[[chr]][(as.numeric(splice_site_coord)) : (as.numeric(splice_site_coord) + intron_padding)])
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))
    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[9]

      fasta_seq <- genome[[chr]][(as.numeric(splice_site_coord) - intron_padding) : (as.numeric(splice_site_coord))]
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))
    }
  } else if (event_type == "dnstream_3ss_exon"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[8]

      fasta_seq <- reverseComplement(genome[[chr]][(as.numeric(splice_site_coord) - exon_padding) : (as.numeric(splice_site_coord))])
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[9]

      fasta_seq <- genome[[chr]][(as.numeric(splice_site_coord)) : (as.numeric(splice_site_coord) + exon_padding)]
      fasta_entry <- data.frame(event = sprintf(">%s:%s", event, event_type),
                                sequence = as.character(fasta_seq))
    }
  }

  return(fasta_entry)
}













##### HOMER ANALYSIS #####
run_directional_HOMER <- function(pos_event, neg_event, output_dir){
  if (!dir.exists(output_dir)) dir.create(output_dir)

  region_list <- c("upstream_5ss_exon",
                   "upstream_5ss_intron",
                   "cassette_3ss_intron",
                   "cassette_3ss_exon",
                   "cassette_5ss_exon",
                   "cassette_5ss_intron",
                   "dnstream_3ss_intron",
                   "dnstream_3ss_exon")
  library(pbmcapply)

  combined_result_df_list <- pbmcapply::pbmclapply(region_list, function(region){
    pos_event_bed <- do.call(rbind, lapply(pos_event, function(x) call_bed_region_from_event(x, region)))
    neg_event_bed <- do.call(rbind, lapply(neg_event, function(x) call_bed_region_from_event(x, region)))

    pos_dir <- sprintf("%s/positive_%s/", output_dir, region)
    neg_dir <- sprintf("%s/negative_%s/", output_dir, region)

    pos_result_list <- run_fg_bg_HOMER(pos_dir, pos_event_bed, neg_event_bed)
    neg_result_list <- run_fg_bg_HOMER(neg_dir, neg_event_bed, pos_event_bed)

    pos_result_df <- filter_HOMER_result(pos_result_list)
    neg_result_df <- filter_HOMER_result(neg_result_list)

    if (!is.null(pos_result_df)){
      pos_result_df$direction <- "positive"
      pos_result_df$region <- region
    }

    if (!is.null(neg_result_df)){
      neg_result_df$direction <- "negative"
      neg_result_df$region <- region
    }
    # annotate


    if (is.null(dim(pos_result_df)) & is.null(dim(neg_result_df))){
      return(NULL)
    } else if (is.null(dim(pos_result_df))) {
      return(neg_result_df)
    } else if (is.null(dim(neg_result_df))) {
      return(pos_result_df)
    } else {
      combined_result_df <- rbind(pos_result_df, neg_result_df)
      return(combined_result_df)
    }
  }, mc.cores = 4)
  all_region_result_df <- do.call(rbind,combined_result_df_list)
  return(all_region_result_df)
}


filter_HOMER_result <- function(result_list){

  filtered_result_list <- list()

  dummy <- lapply(result_list, function(result){


    # check significance
    if (result$Motif_information$fdr < 0.1){
      result_df <- data.frame(consensus = result$Motif_information$consensus,
                              fdr = result$Motif_information$fdr)
      filtered_result_list[[length(filtered_result_list) + 1 ]] <<- result_df
    }
    return()
  })

  combined_filtered_result <- do.call(rbind, filtered_result_list)
  return(combined_filtered_result)
}


run_fg_bg_HOMER <- function(output_dir, event_fg, event_bg, homer_bin_dir = "/Users/harryyang/research/tools/homer/bin", motif_length = 6, num_sim = 100){


  ## write bed file
  # make output dir
  if (!dir.exists(output_dir)) dir.create(output_dir)

  # write bed files
  fg_file <- sprintf("%s/foreground.bed", output_dir)
  bg_file <- sprintf("%s/background.bed", output_dir)

  write.table(event_fg, fg_file, sep = '\t', row.names = F, col.names = F, quote = F)
  write.table(event_bg, bg_file, sep = '\t', row.names = F, col.names = F, quote = F)

  # run perl
  cmd <- sprintf("export PATH=%s:$PATH ; perl %s/findMotifsGenome.pl %s hg19 %s -bg %s -len %s -fdr %s -noknown", homer_bin_dir,  homer_bin_dir, fg_file, output_dir, bg_file, motif_length, num_sim)
  suppressWarnings(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE))

  # read the output
  output_list <- read_denovo_html_results_fdr(output_dir)
  return(output_list)
}


read_denovo_html_results_fdr <- function (path, homer_dir = TRUE)
{
  if (homer_dir == TRUE) {
    path = paste0(path, "/homerResults")
  }
  if (!file.exists(path)) {
    warning(paste("No files found"))
    return(NULL)
  }
  filenames = list.files(path, pattern = "*.info.html")
  df_list = list()
  for (f in filenames) {
    print(f)
    html = readLines(paste(path, f, sep = "/"))
    mypattern = "motif([^<]*).info.html"
    n = gsub(mypattern, "\\1", grep(mypattern, f, value = TRUE))
    df = data.frame(matrix(ncol = 13, nrow = 1))
    colnames(df) = c("motif_name", "consensus", "p_value",
                     "log_p_value", "info_content", "tgt_num", "tgt_pct",
                     "bgd_num", "bgd_pct", "tgt_pos", "bgd_pos", "strand_bias",
                     "multiplicity")
    mypattern = "<H2>Information for ([^<]*)</H2>"
    df$motif_name = gsub(mypattern, "\\1", grep(mypattern,
                                                html, value = TRUE))
    mypattern = paste(".*-([^<]*) \\(Motif ", n, "\\)</H2>",
                      sep = "")
    df$consensus = gsub(mypattern, "\\1", grep(mypattern,
                                               html, value = TRUE))
    mypattern = "<TR><TD>p-value:</TD><TD>([^<]*)</TD></TR>"
    df$p_value = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>log p-value:</TD><TD>([^<]*)</TD></TR>"
    df$log_p_value = gsub(mypattern, "\\1", grep(mypattern,
                                                 html, value = TRUE))
    mypattern = "<TR><TD>log p-value:</TD><TD>([^<]*)</TD></TR>"
    df$log_p_value = gsub(mypattern, "\\1", grep(mypattern,
                                                 html, value = TRUE))
    mypattern = "<TR><TD>FDR:</TD><TD>([^<]*)</TD></TR>"
    df$fdr = gsub(mypattern, "\\1", grep(mypattern,
                                         html, value = TRUE))
    mypattern = "<TR><TD>Information Content per bp:</TD><TD>([^<]*)</TD></TR>"
    df$info_content = gsub(mypattern, "\\1", grep(mypattern,
                                                  html, value = TRUE))
    mypattern = "<TR><TD>Number of Target Sequences with motif</TD><TD>([^<]*)</TD></TR>"
    df$tgt_num = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Percentage of Target Sequences with motif</TD><TD>([^<]*)</TD></TR>"
    df$tgt_pct = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Number of Background Sequences with motif</TD><TD>([^<]*)</TD></TR>"
    df$bgd_num = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Percentage of Background Sequences with motif</TD><TD>([^<]*)</TD></TR>"
    df$bgd_pct = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Average Position of motif in Targets</TD><TD>([^<]*)</TD></TR>"
    df$tgt_pos = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Average Position of motif in Background</TD><TD>([^<]*)</TD></TR>"
    df$bgd_pos = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Strand Bias \\(log2 ratio \\+ to \\- strand density\\)</TD><TD>([^<]*)</TD></TR>"
    df$strand_bias = gsub(mypattern, "\\1", grep(mypattern,
                                                 html, value = TRUE))
    mypattern = "<TR><TD>Multiplicity \\(# of sites on avg that occur together\\)</TD><TD>([^<]*)</TD></TR>"
    df$multiplicity = gsub(mypattern, "\\1", grep(mypattern,
                                                  html, value = TRUE))
    df_list[[n]][["Motif_information"]] = df
    mypattern = "<H4>([^<]*)</H4>"
    length_df = length(gsub(mypattern, "\\1", grep(mypattern,
                                                   html, value = TRUE)))
    df = data.frame(matrix(ncol = 9, nrow = length_df))
    colnames(df) = c("motif_name", "ID", "database", "rank",
                     "score", "offset", "orientation", "original_alignment",
                     "matched_alignment")
    mypattern = "<H4>([^<]*)</H4>"
    df$motif_name = gsub(mypattern, "\\1", grep(mypattern,
                                                html, value = TRUE))
    df <- df %>% tidyr::separate_("motif_name", c("motif_name",
                                                  "ID", "database"), "/", extra = "drop", fill = "right")
    cond <- stringr::str_detect(df$motif_name, "\\(") %>%
      sum(., na.rm = TRUE) > 0
    if (cond == TRUE) {
      df <- df %>% tidyr::separate_("motif_name", c("motif_name",
                                                    "motif_family"), "\\(", extra = "drop", fill = "right")
      df$motif_family <- stringr::str_replace(df$motif_family,
                                              "\\)", "")
    }
    mypattern = "<TR><TD>Match Rank:</TD><TD>([^<]*)</TD></TR>"
    df$rank = gsub(mypattern, "\\1", grep(mypattern, html,
                                          value = TRUE))
    mypattern = "<TR><TD>Score:</TD><TD>([^<]*)</TD</TR>"
    df$score = gsub(mypattern, "\\1", grep(mypattern, html,
                                           value = TRUE))
    mypattern = "<TR><TD>Offset:</TD><TD>([^<]*)</TD</TR>"
    df$offset = gsub(mypattern, "\\1", grep(mypattern, html,
                                            value = TRUE))
    mypattern = "<TR><TD>Orientation:</TD><TD>([^<]*)</TD></TR>"
    df$orientation = gsub(mypattern, "\\1", grep(mypattern,
                                                 html, value = TRUE))
    mypattern = paste(".*-([^<]*) \\(Motif ", n, "\\)</H2>",
                      sep = "")
    df$original_alignment = gsub(mypattern, "\\1", grep(mypattern,
                                                        html, value = TRUE))
    mypattern = ".+>([^<]+)</FONT></TD></TR></TABLE>"
    df$matched_alignment = gsub(mypattern, "\\1", grep(mypattern,
                                                       html, value = TRUE))
    df_list[[n]][["Matches_to_known_motifs"]] = df
  }
  return(df_list)
}



find_motifs_genome_test <- function (x, path, genome, motif_length = c(8, 10, 12), scan_size = 100,
                                     optimize_count = 8, background = "automatic", local_background = FALSE,
                                     only_known = FALSE, only_denovo = FALSE, fdr_num = 0, cores = parallel::detectCores(),
                                     cache = marge:::.calc_free_mem()/4, overwrite = FALSE, keep_minimal = FALSE,
                                     scale_logos = FALSE)
{
  if (overwrite == FALSE & dir.exists(path)) {
    stop("Output directory exists (set `overwrite = TRUE` to bypass)")
  }
  if (background != "automatic" && local_background != FALSE) {
    stop("`background` and `local_background` are mutually exclusive; use only one")
  }
  if (only_known != FALSE & only_denovo != FALSE) {
    stop("Both `only_known` and `only_denovo` set to `TRUE`; pick one")
  }
  if ("data.frame" %in% class(x)) {
    target_bed <- tempfile("target_")
    marge:::.write_bed(x, path = target_bed)
  }
  else {
    if (file.exists(x) != TRUE) {
      stop("Check that your bed file for `x` exists")
    }
    target_bed <- x
  }
  if (!("automatic" %in% background)) {
    if ("data.frame" %in% class(background)) {
      background_bed <- tempfile("background_")
      marge:::.write_bed(background, path = background_bed)
    }
    else {
      if (file.exists(background) != TRUE) {
        stop("Check that your bed file for `background` exists")
      }
      background_bed <- background
    }
  }
  system(paste("mkdir -p", path))
  homer_base <- get_homer_bin()
  cmd <- paste(paste0(homer_base, "findMotifsGenome.pl"),
               target_bed, genome, path, "-len", paste0(motif_length,
                                                        collapse = ","), "-size", scan_size, "-S", optimize_count,
               "-p", cores, "-cache", cache, "-fdr", fdr_num)
  if (!("automatic" %in% background)) {
    cmd <- paste(cmd, "-bg", background_bed)
  }
  if (local_background != FALSE) {
    cmd <- paste(cmd, "-local", local_background)
  }
  if (only_known == TRUE) {
    cmd <- paste(cmd, "-nomotif")
  }
  if (only_denovo == TRUE) {
    cmd <- paste(cmd, "-noknown")
  }
  if (scan_size == "given") {
    cmd <- paste(cmd, "-chopify")
  }
  if (scale_logos == TRUE) {
    cmd <- paste(cmd, "-bits")
  }
  print(cmd)
  system(cmd)
  if (keep_minimal == TRUE) {
    extra_files <- c("homerResults.html", "knownResults.html",
                     "homerMotifs.motifs*", "motifFindingParameters.txt",
                     "seq.autonorm.tsv", "*tmp*")
    extra_dirs <- c("homerResults", "knownResults", "randomizations")
    remove_extra <- paste(c(paste0("rm -f ", path, "/",
                                   extra_files), paste0("rm -Rf ", path, "/", extra_dirs)),
                          collapse = "; ")
    system(remove_extra)
  }
  system("rm -f *.tmp")
}



read_motif_result <- function (path, homer_dir = TRUE)
{
  library(marge)
  if (homer_dir == TRUE) {
    path = paste0(path, "/homerResults")
  }
  if (!file.exists(path)) {
    warning(paste("No files found"))
    return(NULL)
  }
  filenames = list.files(path, pattern = "*.info.html")
  df_list = list()
  for (f in filenames) {
    print(f)
    html = readLines(paste(path, f, sep = "/"))
    mypattern = "motif([^<]*).info.html"
    n = gsub(mypattern, "\\1", grep(mypattern, f, value = TRUE))
    df = data.frame(matrix(ncol = 13, nrow = 1))
    colnames(df) = c("motif_name", "consensus", "p_value",
                     "log_p_value", "info_content", "tgt_num", "tgt_pct",
                     "bgd_num", "bgd_pct", "tgt_pos", "bgd_pos", "strand_bias",
                     "multiplicity")
    mypattern = "<H2>Information for ([^<]*)</H2>"
    df$motif_name = gsub(mypattern, "\\1", grep(mypattern,
                                                html, value = TRUE))
    mypattern = paste(".*-([^<]*) \\(Motif ", n, "\\)</H2>",
                      sep = "")
    df$consensus = gsub(mypattern, "\\1", grep(mypattern,
                                               html, value = TRUE))
    mypattern = "<TR><TD>p-value:</TD><TD>([^<]*)</TD></TR>"
    df$p_value = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>log p-value:</TD><TD>([^<]*)</TD></TR>"
    df$log_p_value = gsub(mypattern, "\\1", grep(mypattern,
                                                 html, value = TRUE))
    mypattern = "<TR><TD>FDR:</TD><TD>([^<]*)</TD></TR>"
    df$FDR = gsub(mypattern, "\\1", grep(mypattern,
                                         html, value = TRUE))
    mypattern = "<TR><TD>Information Content per bp:</TD><TD>([^<]*)</TD></TR>"
    df$info_content = gsub(mypattern, "\\1", grep(mypattern,
                                                  html, value = TRUE))
    mypattern = "<TR><TD>Number of Target Sequences with motif</TD><TD>([^<]*)</TD></TR>"
    df$tgt_num = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Percentage of Target Sequences with motif</TD><TD>([^<]*)</TD></TR>"
    df$tgt_pct = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Number of Background Sequences with motif</TD><TD>([^<]*)</TD></TR>"
    df$bgd_num = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Percentage of Background Sequences with motif</TD><TD>([^<]*)</TD></TR>"
    df$bgd_pct = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Average Position of motif in Targets</TD><TD>([^<]*)</TD></TR>"
    df$tgt_pos = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Average Position of motif in Background</TD><TD>([^<]*)</TD></TR>"
    df$bgd_pos = gsub(mypattern, "\\1", grep(mypattern,
                                             html, value = TRUE))
    mypattern = "<TR><TD>Strand Bias \\(log2 ratio \\+ to \\- strand density\\)</TD><TD>([^<]*)</TD></TR>"
    df$strand_bias = gsub(mypattern, "\\1", grep(mypattern,
                                                 html, value = TRUE))
    mypattern = "<TR><TD>Multiplicity \\(# of sites on avg that occur together\\)</TD><TD>([^<]*)</TD></TR>"
    df$multiplicity = gsub(mypattern, "\\1", grep(mypattern,
                                                  html, value = TRUE))
    df_list[[n]][["Motif_information"]] = df
    mypattern = "<H4>([^<]*)</H4>"
    length_df = length(gsub(mypattern, "\\1", grep(mypattern,
                                                   html, value = TRUE)))
    df = data.frame(matrix(ncol = 9, nrow = length_df))
    colnames(df) = c("motif_name", "ID", "database", "rank",
                     "score", "offset", "orientation", "original_alignment",
                     "matched_alignment")
    mypattern = "<H4>([^<]*)</H4>"
    df$motif_name = gsub(mypattern, "\\1", grep(mypattern,
                                                html, value = TRUE))
    df <- df %>% tidyr::separate_("motif_name", c("motif_name",
                                                  "ID", "database"), "/", extra = "drop", fill = "right")
    cond <- stringr::str_detect(df$motif_name, "\\(") %>%
      sum(., na.rm = TRUE) > 0
    if (cond == TRUE) {
      df <- df %>% tidyr::separate_("motif_name", c("motif_name",
                                                    "motif_family"), "\\(", extra = "drop", fill = "right")
      df$motif_family <- stringr::str_replace(df$motif_family,
                                              "\\)", "")
    }
    mypattern = "<TR><TD>Match Rank:</TD><TD>([^<]*)</TD></TR>"
    df$rank = gsub(mypattern, "\\1", grep(mypattern, html,
                                          value = TRUE))
    mypattern = "<TR><TD>Score:</TD><TD>([^<]*)</TD</TR>"
    df$score = gsub(mypattern, "\\1", grep(mypattern, html,
                                           value = TRUE))
    mypattern = "<TR><TD>Offset:</TD><TD>([^<]*)</TD</TR>"
    df$offset = gsub(mypattern, "\\1", grep(mypattern, html,
                                            value = TRUE))
    mypattern = "<TR><TD>Orientation:</TD><TD>([^<]*)</TD></TR>"
    df$orientation = gsub(mypattern, "\\1", grep(mypattern,
                                                 html, value = TRUE))
    mypattern = paste(".*-([^<]*) \\(Motif ", n, "\\)</H2>",
                      sep = "")
    df$original_alignment = gsub(mypattern, "\\1", grep(mypattern,
                                                        html, value = TRUE))
    mypattern = ".+>([^<]+)</FONT></TD></TR></TABLE>"
    df$matched_alignment = gsub(mypattern, "\\1", grep(mypattern,
                                                       html, value = TRUE))
    df_list[[n]][["Matches_to_known_motifs"]] = df
  }
  return(df_list)
}



call_bed_region_from_event <- function(event, event_type, intron_padding = 300, exon_padding = 50){
  bed_region_list <- strsplit(event, ":")[[1]]
  chr <- bed_region_list[3]
  direction <- bed_region_list[4]
  if (event_type == "upstream_5ss_exon"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[9]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord),
                              end = as.numeric(splice_site_coord) + exon_padding)

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[8]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord) - exon_padding,
                              end = as.numeric(splice_site_coord))

    }
  } else if (event_type == "upstream_5ss_intron"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[9]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord) - intron_padding,
                              end = as.numeric(splice_site_coord))

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[8]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord),
                              end = as.numeric(splice_site_coord) + intron_padding)

    }
  } else if (event_type == "cassette_3ss_intron"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[6]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord),
                              end = as.numeric(splice_site_coord) + intron_padding)

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[5]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord) - intron_padding,
                              end = as.numeric(splice_site_coord))

    }
  } else if (event_type == "cassette_3ss_exon"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[6]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord) - exon_padding,
                              end = as.numeric(splice_site_coord))

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[5]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord),
                              end = as.numeric(splice_site_coord) + exon_padding)

    }
  } else if (event_type == "cassette_5ss_exon"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[5]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord),
                              end = as.numeric(splice_site_coord) + exon_padding)

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[6]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord) - exon_padding,
                              end = as.numeric(splice_site_coord))

    }
  } else if (event_type == "cassette_5ss_intron"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[5]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord) - intron_padding,
                              end = as.numeric(splice_site_coord))
    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[6]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord),
                              end = as.numeric(splice_site_coord) + intron_padding)
    }
  } else if (event_type == "dnstream_3ss_intron"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[8]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord),
                              end = as.numeric(splice_site_coord) + intron_padding)

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[9]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord) - intron_padding,
                              end = as.numeric(splice_site_coord))

    }
  } else if (event_type == "dnstream_3ss_exon"){
    if (direction == "-"){
      splice_site_coord <- bed_region_list[8]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord) - exon_padding,
                              end = as.numeric(splice_site_coord))

    } else if (direction == "+"){
      splice_site_coord <- bed_region_list[9]

      bed_entry <- data.frame(chr = chr,
                              start = as.numeric(splice_site_coord),
                              end = as.numeric(splice_site_coord) + exon_padding)

    }
  }
  bed_entry$id <- event
  bed_entry$category <- event_type
  bed_entry$strand <- direction

  return(bed_entry)
}

