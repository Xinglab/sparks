
kd_library <- readRDS("/Users/harryyang/transfer/Clean_minfiltered.SE.library.rds")

##### MATSUMOTO #####
matsumoto_sparks_file <- '//Users/harryyang/transfer/Spl_inhibitor/Matsumoto_PC3_Pladienolide.SPARKS.rds'
matsumoto_sparks <- readRDS(matsumoto_sparks_file)
# matsumoto_mats <- import_SPARKS_MATS_for_analysis(matsumoto_sparks, "SE")

matsumoto_result <- matsumoto_sparks@SPARKS_analysis_result$SE

# generate bar plot
add_plot_title(generate_enrichment_barplot(matsumoto_result,
                                           bar_color = "goldenrod2",
                                           num_plot = 10),
               "PC3 Pladienolide-B",
               title_color = "goldenrod2") # 350 x 400

##### POGACAR #####
pogacar_sparks_file <- '/Users/harryyang/transfer/Spl_inhibitor/Pogacar_A549_Indisulam.SPARKS.rds'
pogacar_sparks <- readRDS(pogacar_sparks_file)
# pogacar_mats <- import_SPARKS_MATS_for_analysis(pogacar_sparks, "SE")

pogacar_result <- pogacar_sparks@SPARKS_analysis_result$SE

# generate bar plot
add_plot_title(generate_enrichment_barplot(pogacar_result,
                                           bar_color = "sandybrown",
                                           num_plot = 10),
               "A549 Indisulam",
               title_color = "sandybrown") # 350 x 400

##### WU #####
wu_sparks_file <- '/Users/harryyang/transfer/Spl_inhibitor/Rh18_Sudemycine_D1.SPARKS.rds'
wu_sparks <- readRDS(wu_sparks_file)
# wu_mats <- import_SPARKS_MATS_for_analysis(wu_sparks, "SE")

wu_result <- wu_sparks@SPARKS_analysis_result$SE

# generate bar plot
add_plot_title(generate_enrichment_barplot(wu_result,
                                           bar_color = "dodgerblue1",
                                           num_plot = 10),
               "LNCaP ESRP1/2 siRNA",
               title_color = "dodgerblue2") # 350 x 400
