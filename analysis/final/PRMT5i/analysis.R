##### INITIALIZATION #####
library(SPARKS)
library(ggplot2)
library(dplyr)



kd_library <- readRDS("/Users/harryyang/transfer/Clean_minfiltered.SE.library.rds")

##### ANALYSIS #####
# read in SNRPB data
snrpb_sparks_file <- '/Users/harryyang/transfer/PRMT5i/Correa_SNRPB_siRNA.SPARKS.rds'
snrpb_sparks <- readRDS(snrpb_sparks_file)
snrpb_mats <- import_SPARKS_MATS_for_analysis(snrpb_sparks, "SE")

snrpb_result <- snrpb_sparks@SPARKS_analysis_result$SE

add_plot_title(generate_enrichment_barplot(snrpb_result,
                                           bar_color = "lavender",
                                           num_plot = 10),
               "U251 SNRPB siRNA",
               title_color = "mediumpurple2") # 350 x 400


# process braun data
braun_sparks_file <- '/Users/harryyang/transfer/PRMT5i/Braun_PRMT5i_2017.SPARKS.rds'
braun_sparks <- readRDS(braun_sparks_file)
braun_mats <- import_SPARKS_MATS_for_analysis(braun_sparks, "SE")

braun_result <- braun_sparks@SPARKS_analysis_result$SE

# add SNRPB data
braun_new_result <- add_custom_library_to_SPARKS_test_result(braun_mats,
                                                             braun_result,
                                                             snrpb_mats,
                                                             "U251_SNRPB_siRNA",
                                                             kd_library)
# generate bar plot
add_plot_title(generate_enrichment_barplot(braun_new_result,
                                           bar_color = "mediumpurple2",
                                           num_plot = 10),
               "U87 PRMT5 inhibition",
               title_color = "mediumpurple2") # 350 x 400


# process a427 data
a427_sparks_file <- '/Users/harryyang/transfer/PRMT5i/A427_PRMT5_inhibition.30.SPARKS.rds'
a427_sparks <- readRDS(a427_sparks_file)
a427_mats <- import_SPARKS_MATS_for_analysis(a427_sparks, "SE")

a427_result <- a427_sparks@SPARKS_analysis_result$SE

# add SNRPB data
a427_new_result <- add_custom_library_to_SPARKS_test_result(a427_mats,
                                                            a427_result,
                                                            snrpb_mats,
                                                            "U251_SNRPB_siRNA",
                                                            kd_library)
# generate bar plot
add_plot_title(generate_enrichment_barplot(a427_new_result,
                                           bar_color = "mediumpurple2",
                                           num_plot = 10),
               "A427 PRMT5 inhibition",
               title_color = "mediumpurple2") # 350 x 400


# process h1975 data
h1975_sparks_file <- '/Users/harryyang/transfer/PRMT5i/H1975_PRMT5_inhibition.30.SPARKS.rds'
h1975_sparks <- readRDS(h1975_sparks_file)
h1975_mats <- import_SPARKS_MATS_for_analysis(h1975_sparks, "SE")

h1975_result <- h1975_sparks@SPARKS_analysis_result$SE

# add SNRPB data
h1975_new_result <- add_custom_library_to_SPARKS_test_result(h1975_mats,
                                                             h1975_result,
                                                             snrpb_mats,
                                                             "U251_SNRPB_siRNA",
                                                             kd_library)
# generate bar plot
add_plot_title(generate_enrichment_barplot(h1975_new_result,
                                           bar_color = "mediumpurple2",
                                           num_plot = 10),
               "H1975 PRMT5 inhibition",
               title_color = "mediumpurple2") # 350 x 400

# process h441 data
h441_sparks_file <- '/Users/harryyang/transfer/PRMT5i/H441_PRMT5_inhibition.30.SPARKS.rds'
h441_sparks <- readRDS(h441_sparks_file)
h441_mats <- import_SPARKS_MATS_for_analysis(h441_sparks, "SE")

h441_result <- h441_sparks@SPARKS_analysis_result$SE

# add SNRPB data
h441_new_result <- add_custom_library_to_SPARKS_test_result(h441_mats,
                                                            h441_result,
                                                            snrpb_mats,
                                                            "U251_SNRPB_siRNA",
                                                            kd_library)
# generate bar plot
add_plot_title(generate_enrichment_barplot(h441_new_result,
                                           bar_color = "mediumpurple2",
                                           num_plot = 10),
               "H441 PRMT5 inhibition",
               title_color = "mediumpurple2") # 350 x 400

# process k562 data
k562_sparks_file <- '/Users/harryyang/transfer/PRMT5i/K562_PRMT5_inhibition.30.SPARKS.rds'
k562_sparks <- readRDS(k562_sparks_file)
k562_mats <- import_SPARKS_MATS_for_analysis(k562_sparks, "SE")

k562_result <- k562_sparks@SPARKS_analysis_result$SE

# add SNRPB data
k562_new_result <- add_custom_library_to_SPARKS_test_result(k562_mats,
                                                            k562_result,
                                                            snrpb_mats,
                                                            "U251_SNRPB_siRNA",
                                                            kd_library)
# generate bar plot
add_plot_title(generate_enrichment_barplot(k562_new_result,
                                           bar_color = "mediumpurple2",
                                           num_plot = 10),
               "K562 PRMT5 inhibition",
               title_color = "mediumpurple2") # 350 x 400


# process maron data
maron_sparks_file <- '/Users/harryyang/transfer/PRMT5i/Maron_PRMT5i_4days.SPARKS.rds'
maron_sparks <- readRDS(maron_sparks_file)
maron_mats <- import_SPARKS_MATS_for_analysis(maron_sparks, "SE")

maron_result <- maron_sparks@SPARKS_analysis_result$SE

# add SNRPB data
maron_new_result <- add_custom_library_to_SPARKS_test_result(maron_mats,
                                                             maron_result,
                                                             snrpb_mats,
                                                             "U251_SNRPB_siRNA",
                                                             kd_library)
# generate bar plot
add_plot_title(generate_enrichment_barplot(maron_new_result,
                                           bar_color = "mediumpurple2",
                                           num_plot = 10),
               "A549 PRMT5 inhibition",
               title_color = "mediumpurple2") # 350 x 400


# process rad data
rad_sparks_file <- '/Users/harryyang/transfer/PRMT5i/Rad_PRMT5i_2019.SPARKS.rds'
rad_sparks <- readRDS(rad_sparks_file)
rad_mats <- import_SPARKS_MATS_for_analysis(rad_sparks, "SE")

rad_result <- rad_sparks@SPARKS_analysis_result$SE

# add SNRPB data
rad_new_result <- add_custom_library_to_SPARKS_test_result(rad_mats,
                                                           rad_result,
                                                           snrpb_mats,
                                                           "U251_SNRPB_siRNA",
                                                           kd_library)
# generate bar plot
add_plot_title(generate_enrichment_barplot(rad_new_result,
                                           bar_color = "mediumpurple2",
                                           num_plot = 10),
               "THP-1 PRMT5 inhibition",
               title_color = "mediumpurple2") # 350 x 400


