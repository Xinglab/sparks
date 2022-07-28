library('SPARKS')
library("optparse")
library(maftools)
library(ggplot2)
library(stringr)
library(reshape2)
library(RColorBrewer)
# load gtools for mixedsort for chromosome-level coordinate sorting
library(gtools)
library(data.table)
library(R.utils)
library(dplyr)
library(pbmcapply)

option_list = list(
  make_option(c("--cancer_type"), type = "character", default = NULL, help = "Cancer Type"),

  make_option(c("--psi_SE"), type = "character", default = NULL, help = "PSI matrix for SE"),
  make_option(c("--psi_A3SS"), type = "character", default = NULL, help = "PSI matrix for A3SS"),
  make_option(c("--psi_A5SS"), type = "character", default = NULL, help = "PSI matrix for A5SS"),
  make_option(c("--psi_RI"), type = "character", default = NULL, help = "PSI matrix for RI"),

  make_option(c("--exp"), type = "character", default = NULL, help = "Expression matrix"),

  make_option(c("--CLIP_SE"), type = "character", default = NULL, help = "ENCODE CLIP-Seq overlap for SE events"),
  make_option(c("--CLIP_A3SS"), type = "character", default = NULL, help = "ENCODE CLIP-Seq overlap for A3SS events"),
  make_option(c("--CLIP_A5SS"), type = "character", default = NULL, help = "ENCODE CLIP-Seq overlap for A5SS events"),
  make_option(c("--CLIP_RI"), type = "character", default = NULL, help = "ENCODE CLIP-Seq overlap for RI events"),

  make_option(c("--exon_anno_SE"), type = "character", default = NULL, help = "Exon annotation for SE events"),
  make_option(c("--exon_anno_A3SS"), type = "character", default = NULL, help = "Exon annotation for A3SS events"),
  make_option(c("--exon_anno_A5SS"), type = "character", default = NULL, help = "Exon annotation for A5SS events"),
  make_option(c("--exon_anno_RI"), type = "character", default = NULL, help = "Exon annotation for RI events"),

  make_option(c("--MATS_SE"), type = "character", default = NULL, help = "MATS result for SE events"),
  make_option(c("--MATS_A3SS"), type = "character", default = NULL, help = "MATS result for A3SS events"),
  make_option(c("--MATS_A5SS"), type = "character", default = NULL, help = "MATS result for A5SS events"),
  make_option(c("--MATS_RI"), type = "character", default = NULL, help = "MATS result for RI events"),
  
  make_option(c("--SPARKS_library"), type = "character", default = NULL, help = "Library for SPARKS")
)

# TODO - add options to run SPARKS analysis by default

opt <- parse_args(OptionParser(option_list = option_list))

input_psi_file_SE <- opt$psi_SE
input_psi_file_A3SS <- opt$psi_A3SS
input_psi_file_A5SS <- opt$psi_A5SS
input_psi_file_RI <- opt$psi_RI

sparks_obj <- new("SPARKS", study = opt$cancer_type)
sparks_obj <- import_PSI_df(sparks_obj, input_psi_file_SE, "SE")
sparks_obj <- import_PSI_df(sparks_obj, input_psi_file_A3SS, "A3SS")
sparks_obj <- import_PSI_df(sparks_obj, input_psi_file_A5SS, "A5SS")
sparks_obj <- import_PSI_df(sparks_obj, input_psi_file_RI, "RI")


# import expression df
if (!(is.null(opt$exp))){
  input_exp_file <- opt$exp
  sparks_obj <- import_expression_matrix(sparks_obj, input_exp_file)
}


# import clip data
if (!(is.null(opt$CLIP_SE))){
  input_clip_SE <- opt$CLIP_SE
  sparks_obj <- import_ENCODE_CLIP_intersect_data(sparks_obj, input_clip_SE, "SE")
}

if (!(is.null(opt$CLIP_A3SS))){
  input_clip_A3SS <- opt$CLIP_A3SS
  sparks_obj <- import_ENCODE_CLIP_intersect_data(sparks_obj, input_clip_A3SS, "A3SS")
}

if (!(is.null(opt$CLIP_A5SS))){
  input_clip_A5SS <- opt$CLIP_A5SS
  sparks_obj <- import_ENCODE_CLIP_intersect_data(sparks_obj, input_clip_A5SS, "A5SS")
}

if (!(is.null(opt$CLIP_RI))){
  input_clip_RI <- opt$CLIP_RI
  sparks_obj <- import_ENCODE_CLIP_intersect_data(sparks_obj, input_clip_RI, "RI")
}

# import exon annotation data for downstream analysis
input_exon_anno_SE <- opt$exon_anno_SE
sparks_obj <- import_exon_annotation(sparks_obj, input_exon_anno_SE, 'SE')
input_exon_anno_A3SS <- opt$exon_anno_A3SS
sparks_obj <- import_exon_annotation(sparks_obj, input_exon_anno_A3SS, 'A3SS')
input_exon_anno_A5SS <- opt$exon_anno_A5SS
sparks_obj <- import_exon_annotation(sparks_obj, input_exon_anno_A5SS, 'A5SS')
input_exon_anno_RI <- opt$exon_anno_RI
sparks_obj <- import_exon_annotation(sparks_obj, input_exon_anno_RI, 'RI')

# MOTIF analysis
if (!(is.null(opt$motif_SE))){
  input_motif_dir_SE <- opt$motif_SE
  sparks_obj <- import_motif_count_data(sparks_obj, input_motif_dir_SE, 'SE')
}


# import MATS result
if (!(is.null(opt$MATS_SE))){
  input_MATS_SE <- opt$MATS_SE
  sparks_obj <- import_MATS(sparks_obj, input_MATS_SE, 'SE')
}
if (!(is.null(opt$MATS_A3SS))){
  input_MATS_A3SS <- opt$MATS_A3SS
  sparks_obj <- import_MATS(sparks_obj, input_MATS_A3SS, 'A3SS')
}
if (!(is.null(opt$MATS_A5SS))){
  input_MATS_A5SS <- opt$MATS_A5SS
  sparks_obj <- import_MATS(sparks_obj, input_MATS_A5SS, 'A5SS')
}
if (!(is.null(opt$MATS_RI))){
  input_MATS_RI <- opt$MATS_RI
  sparks_obj <- import_MATS(sparks_obj, input_MATS_RI, 'RI')
}

# perform SPARKS analysis if the option is used
if (!(is.null(opt$SPARKS_library))){
  kd_library_all <- readRDS(opt$SPARKS_library)
  analysis_result <- perform_SPARKS_analysis_for_all_splice_types(sparks_obj, 
    kd_library_all, 
    test_study = opt$cancer_type)
  print("Performing SPARKS analysis")
  sparks_obj$SPARKS_analysis_result <<- analysis_result
}


# save the object
saveRDS(sparks_obj, sprintf("%s.SPARKS.rds", opt$cancer_type))




