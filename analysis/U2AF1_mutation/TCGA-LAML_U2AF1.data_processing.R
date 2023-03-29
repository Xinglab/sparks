##### CODE #####
library(SPARKS)
library(ggplot2)
score_method <- "GSEA"
spl_types <- c("SE")



##### MUTATION ANALYSIS ######
library(maftools)

laml_maf <- maftools::read.maf('/Users/harryyang/research/Xing_Lab/sparks/data/SF_mutation/MC3.TCGA-LAML.SNP.maf.gz',
                               isTCGA = T)

# geneate lollipop plot
lollipopPlot(laml_maf, gene = "U2AF1", AACol = "HGVSp", labelPos = "all")

# extract tumor samples with WXS data
laml_wxs_sample_list <- as.character(laml_maf@variants.per.sample$Tumor_Sample_Barcode)

# extract samples with s34f mutation
laml_u2af1_subset <- subsetMaf(laml_maf, genes = "U2AF1")
laml_wxs_s34f_sample_list <- as.character(laml_u2af1_subset@data[laml_u2af1_subset@data$HGVSp_Short == "p.S34F", ]$Tumor_Sample_Barcode)

# extract samples with WT U2AF1
laml_wxs_mut_sample_list <- as.character(laml_u2af1_subset@data$Tumor_Sample_Barcode)
laml_wxs_wt_sample_list <- laml_wxs_sample_list[!(laml_wxs_sample_list %in% laml_wxs_mut_sample_list)]


## read sample list - this is the order from the mats
laml_sample_list <- colnames(read.csv('/Users/harryyang/research/Xing_Lab/sparks/data/SF_mutation/TCGA-LAML/TCGA_LAML_RL50_rmatspost_list.txt', sep = ',', check.names = F))

# note - the last part is trimmed as it might not match between WXS and RNA-Seq,
# so tumor samples are mached by the first 10
laml_sample_list_clean <- unlist(lapply(laml_sample_list, function(x) paste(strsplit(strsplit(strsplit(x, "/")[[1]][8], "\\.")[[1]][1], "-")[[1]][1:3], collapse = "-")))

# laml_s34_sample_list <- c("TCGA-AB-2996-03A",
#                           "TCGA-AB-2847-03A",
#                           "TCGA-AB-2843-03A",
#                           "TCGA-AB-2912-03A")
#
# # these two are S34Y
# laml_s34y_sample_list <- c( "TCGA-AB-2882-03A",
#                             "TCGA-AB-2861-03A")
#
# laml_s34_mut_index <- which(laml_sample_list_clean %in% laml_s34_sample_list)
# laml_s34_wt_index <- which(!(laml_sample_list_clean %in% union(laml_s34_sample_list, laml_s34y_sample_list)))

laml_s34_mut_index <- which(laml_sample_list_clean %in% laml_wxs_s34f_sample_list)
laml_s34_wt_index <- which(laml_sample_list_clean %in% laml_wxs_wt_sample_list)


laml_study <- "TCGA-LAML_U2AF1_S34F"

### Run SPARKS analysis
spl_type = "SE"
# filter raw splicing data
input_mats_file <- sprintf('/Users/harryyang/research/Xing_Lab/sparks/data/SF_mutation/TCGA-LAML//TCGA-LAML.%s.MATS_df.txt', spl_type)
input_mats_df <- as.data.frame(data.table::fread(input_mats_file))

# stratify the samples by the mutation status
study_mats <- calculate_new_psi_for_subset_of_samples(input_mats_df, laml_s34_wt_index, laml_s34_mut_index, min_threshold = 10) # KD parity check

# cleanup the names to match the library
study_mats$event <- unlist(lapply(study_mats$event, function(x) rewrite_event_coordinates(x)))

# sort by beta and remove NA, possibly generated from the subset operation
study_mats <- study_mats %>% na.omit() %>% arrange(-beta)

# save the result for downstream stuffs
data.table::fwrite(study_mats, "~/transfer/SF_mutation/TCGA-LAML.SE.U2AF1_S34F_sorted.MATS.df.txt")



# make table for sashimi
bam_list <- data.table::fread("~/test_17", header = F)

laml_s34_mut_list <- laml_sample_list_clean[(laml_sample_list_clean %in% laml_wxs_s34f_sample_list)]
laml_s34_wt_list <- laml_sample_list_clean[(laml_sample_list_clean %in% laml_wxs_wt_sample_list)]

mut_bam_list <- unlist(lapply(laml_s34_mut_list, function(sample) bam_list$V1[grep(sample, bam_list$V1)]))
wt_bam_list <- unlist(lapply(laml_s34_wt_list, function(sample) bam_list$V1[grep(sample, bam_list$V1)]))

# write table
sashimi_table <- data.frame(sample = c(mut_bam_list,
                                       wt_bam_list),
                            group = c(rep("U2AF1 S34F", length(mut_bam_list)),
                                      rep("U2AF1 WT", length(wt_bam_list))))
data.table::fwrite(sashimi_table, "~/sample_table.txt")
