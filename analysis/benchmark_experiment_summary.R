##### CODE #####
library(SPARKS)
library(ggplot2)
spl_types <- c("SE", "A3SS", "A5SS")

##### ANALYSIS #####
merged_kd_library <- readRDS('/Users/harryyang/research/Xing_Lab/sparks/data/ENCODE_library.KD_and_KO.count_20.rds')

# determine replicate RBP experiments
shrna_study_list <- names(merged_kd_library$SE)[grep("shRNA", names(merged_kd_library$SE))]

# generate study -> RBP table
shrna_rbp_info_df <- data.frame(study = shrna_study_list,
                                rbp = unlist(lapply(shrna_study_list, function(x) strsplit(x, "_")[[1]][2])))
shrna_rbp_num_exp <- table(shrna_rbp_info_df$rbp)
shrna_rbp_kd2_list <- names(shrna_rbp_num_exp)[shrna_rbp_num_exp == 2]

# determine replicate RBP experiments
crispr_study_list <- names(merged_kd_library$SE)[grep("CRISPR", names(merged_kd_library$SE))]

# generate study -> RBP table
crispr_rbp_info_df <- data.frame(study = crispr_study_list,
                                 rbp = unlist(lapply(crispr_study_list, function(x) strsplit(x, "_")[[1]][2])))
crispr_rbp_num_exp <- table(crispr_rbp_info_df$rbp)
crispr_rbp_kd2_list <- names(crispr_rbp_num_exp)[crispr_rbp_num_exp == 2]

##### SUMMARY TABLE #####
num_shrna_k562 <- length(grep(shRNA_study_list, pattern = "K562"))
num_shrna_hepg2 <- length(grep(shRNA_study_list, pattern = "HepG2"))
num_crispr_k562 <- length(grep(shRNA_study_list, pattern = "K562"))
num_crispr_hepg2 <- length(grep(shRNA_study_list, pattern = "HepG2"))

num_shrna_common <- length(shrna_rbp_kd2_list)
num_crispr_common <- length(crispr_rbp_kd2_list)

num_rbp_in_shrna_crispr <- length(intersect(unique(shrna_rbp_info_df$rbp), unique(crispr_rbp_info_df$rbp)))


