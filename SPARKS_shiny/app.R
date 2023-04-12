#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize=1200*1024^2) # set max size at 1.2g
library(shiny)
library(shinyWidgets)
library(dplyr)
library(SPARKS)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("SPARKS"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            tabsetPanel(type = "tabs",
                        id = "sidetabs",
                        tabPanel(title = "Plot",
                                 value = "tab_plot",
                                 # Copy the line below to make a file upload manager
                                 fileInput("sparks_rds",
                                           label = h5("Processed SPARKS Object"),
                                           accept = ".rds"),
                                 uiOutput("plot_button"),
                                 uiOutput("numplot_slider"),
                                 uiOutput("custom_highlight_selection"),
                                 uiOutput("custom_library_button")),
                        tabPanel(title = "Custom Signatures",
                                 value = "tab_custom_signature",
                                 # uiOutput("ref_library_file"),
                                 uiOutput("custom_library_file"),
                                 uiOutput("custom_library_process_button"),
                                 uiOutput("reset_custom_library_process"))  # reset the process
            )

        ),

        # Show a plot of the generated distribution
        mainPanel(
            h4("Enrichment Barplot"),
            plotOutput("barplot", height = '600px')
        )
    ),
    mainPanel(h4("SPARKS Analysis Result Table"),
              dataTableOutput("sparks_table"),
              width = 12)
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    # define necessary data variables
    input_mats <- reactiveVal(NULL)
    input_sparks_result <- reactiveVal(NULL)
    # reference_library <- reactiveVal(NULL)

    # add button to generate plot
    output$plot_button <- renderUI({
        req(input$sparks_rds)
        # req(input_sparks)
        actionButton(inputId = "plot_gen",
                     label = "Generate SPARKS barplot and table")
    })

    # add num plot slider
    output$numplot_slider <- renderUI({
        req(input$plot_gen)
        tagList(hr(),
                sliderInput(inputId = 'num_plot',
                            label = '# Experiments to plot',
                            min = 5, max = 50, value = 10)
        )
    })

    # # add custom highlight button
    # output$custom_highlight_button <- renderUI({
    #     req(input$plot_gen)
    #     actionButton(inputId = "custom_highlight_activated",
    #                  label = "Highlight specific RBPs")
    # })

    # select specific RBPs to highlight
    output$custom_highlight_selection <- renderUI({
        req(input$plot_gen)

        # query rbp list in the dataset
        unique_rbp_list <- sort(unique(unlist(lapply(input_sparks_result()$S1,
                                                     function(x) strsplit(x, "_")[[1]][2]))))

        # make options
        tagList(hr(),
                pickerInput("custom_highlight_rbps",
                            "RBPs to highlight",
                            choices = unique_rbp_list,
                            options = list(`actions-box` = TRUE),
                            multiple = T)
        )
    })


    # BELOW LINE IS FOR FULL LIBRARY LOAD - not necessary for now
    #     # select library file
    #     output$ref_library_file <- renderUI({
    #         req(input$custom_library_activated)
    #         print("AA")
    #         tagList(fileInput(inputId = "ref_library_rds",
    #                           label = "Reference Library",
    #                           accept = ".rds"),
    #                 actionButton(inputId = "ref_library_load",
    #                              label = "Load the reference library"))
    #     })
    #
    #     # load library file
    #     ref_library_loaded <- observeEvent(input$ref_library_load, {
    #
    #         print("BB")
    #         withProgress(message = 'Importing reference library data', value = 0, {
    #             print(input$ref_library_rds)
    #             # read the rds
    #             ref_library <- readRDS(input$ref_library_rds$datapath)
    #             incProgress(0.5, detail = "Imported data")
    #
    #             # update the analysis result with new result
    #             reference_library(ref_library)
    #             incProgress(0.5, detail = "Updated the system")
    #
    #         })
    #         return(TRUE)
    #     })

    ##### CUSTOM LIBRARY #####
    # add custom highlight button
    output$custom_library_button <- renderUI({
        req(input$plot_gen)
        actionButton(inputId = "custom_library_activated",
                     label = "Add custom library into analysis")
    })

    # move to the custom signature tab if the button is clicked
    observeEvent(input$custom_library_activated, {
        updateTabsetPanel(session, "sidetabs",
                          selected = 'tab_custom_signature')
    })

    # display widgets for the custom signatures
    output$custom_library_file <- renderUI({
        print("XX")
        # if(!req(input$custom_library_activated, cancelOutput = T)){  # if the custom library is not activated

        if(is.null(input$custom_library_activated)){  # if the custom library is not activated
            print("AAB")
            h5("This will be available after loading the analysis data")

        } else {  # show input fields for custom library analysis
            tagList(fileInput(inputId = "custom_library_rds",
                              label = "Custom Signature",
                              accept = ".rds",
                              placeholder = "SPARKS object"),
                    textInput(inputId = "custom_library_study",
                              label = "Unique signature name",
                              placeholder = "CellLine_RBP_perturbation"))
        }

    })

    # generate process button when the necessary parts are ready
    output$custom_library_process_button <- renderUI({
        req(input$custom_library_rds, input$custom_library_study)
        actionButton(inputId = "custom_library_process_activated",
                     label = "Process custom library")

    })

    # store last imported library to reduce possibility of errors importing the same file
    # - this can happen when user clicks import twice so just skip it all together
    last_imported_custom_file_temploc <- reactiveVal(NULL)

    # store the status of the custom library process for repeated run
    custom_process_finished <- reactiveVal(FALSE)

    # read and run the analysis using this specific library
    custom_mats <- observeEvent(input$custom_library_process_activated, ignoreInit = TRUE, {

        # set the status to false
        custom_process_finished(FALSE)

        # locate the temp file
        custom_lib_path <- input$custom_library_rds$datapath

        # check if it matches the last imported file - meaning it shouldn't be imported
        if (!is.null(last_imported_custom_file_temploc())){
            if (as.character(custom_lib_path) == as.character(last_imported_custom_file_temploc())){  # if matches
                # note the temploc can be NULL, so using '==' to compare doesn't work if not null checked
                # change the button to signify this
                updateActionButton(inputId = "custom_library_process_activated",
                                   label = "Signature already imported")

                # use validate to stop the run - TODO - this should be improved
                validate(
                    need(custom_lib_path != last_imported_custom_file_temploc(),
                         message = "Signature already imported")
                )
            }
        } # if it doesn't match, we will update the temploc after the validation ends

        # check if the signature name is already in the data
        # - this will throw error in downstream processing and confuse users
        if (input$custom_library_study %in% unique(input_sparks_result()$S1)){  # if it's already used
            # change the button to signify this
            updateActionButton(inputId = "custom_library_process_activated",
                               label = "Duplicate signature name")

            # use validate to stop the run - TODO - this should be improved
            validate(
                need(!(input$custom_library_study %in% unique(input_sparks_result()$S1)),
                     message = "Please choose a unique study name")
            )
        }

        # update the button
        updateActionButton(inputId = "custom_library_process_activated",
                           label = "Processing...")

        # import with the progress
        withProgress(message = 'Processing custom library data', value = 0, {
            # read the rds
            custom_library_sparks <- readRDS(custom_lib_path)
            incProgress(0.2, detail = "Imported data")
            print("A")
            # get the mats for custom study
            custom_library_mats <- import_SPARKS_MATS_for_analysis(custom_library_sparks,
                                                                   "SE")
            incProgress(0.2, detail = "Imported MATS for analysis")
            print("B")

            # import the count data
            # - rather than needing to load the whole library, we only use count info
            # - as this is all the necessary info from the library for the new run
            event_count_df <- data.table::fread("./SPARKS_shiny/data/event_count.SE.df.txt")
            incProgress(0.2, detail = "Imported event count df for analysis")


            # run custom analysis using this
            new_result <- add_custom_library_to_SPARKS_test_result(input_mats(),
                                                                   input_sparks_result(),
                                                                   custom_library_mats,
                                                                   input$custom_library_study,
                                                                   event_count_df = event_count_df)
            incProgress(0.2, detail = "Performed SPARKS analysis")
            print("C")

            # update the analysis result with new result
            input_sparks_result(new_result)
            incProgress(0.2, detail = "Updated analysis result")
            print("D")

        })

        # update the last import
        last_imported_custom_file_temploc(custom_lib_path)

        # update the button
        updateActionButton(inputId = "custom_library_process_activated", label = "Process finished")

        # change the status to complete to make reset button
        custom_process_finished(TRUE)
    })

    # make button to reset this process if the previous process is finished
    output$reset_custom_library_process <- renderUI({
        req(custom_process_finished() == TRUE)
        actionButton(inputId = "reset_custom_library_process",
                     label = "Import another custom signature")
    })

    # if the reset button is pressed, make the UI again
    observeEvent(input$reset_custom_library_process, {
        # update the flag to FALSE, as new dataset is imported in the new slot
        custom_process_finished(FALSE)
        # make the UI again
        print("XXZZ")
        output$custom_library_file <- renderUI({

            tagList(fileInput(inputId = "custom_library_rds",
                              label = "Custom Signature",
                              accept = ".rds",
                              placeholder = "SPARKS object"),
                    textInput(inputId = "custom_library_study",
                              label = "Unique signature name",
                              placeholder = "CellLine_RBP_perturbation"))
        })
    })



    ##### CORE FUNCTIONALITY #####
    # Import the input SPARKS object
    input_sparks <- observeEvent(input$plot_gen, ignoreInit = TRUE, {
        withProgress(message = 'Importing data', value = 0, {
            print(input$sparks_rds$datapath)
            input_sparks_obj <- readRDS(input$sparks_rds$datapath)
            incProgress(0.3, detail = "Imported data")

            # import necessary information
            input_sparks_mats <-import_SPARKS_MATS_for_analysis(input_sparks_obj, 'SE')
            input_result <- input_sparks_obj@SPARKS_analysis_result$SE
            incProgress(0.3, detail = "Processed data")

            # update the variables
            input_mats(input_sparks_mats)
            input_sparks_result(input_result)
            incProgress(0.3, detail = "Updated the necessary values")

        })
        return(input_sparks_obj)
    })


    # show table
    output$sparks_table <- renderDataTable({
        req(input$plot_gen)
        result <- input_sparks_result()

        # calculate FWER
        result_padj <- calculate_SPARKS_padj(result)

        # only keep relevant information for showing
        relevant_cols <- c("S1",
                           "gsea_pos_score",
                           "gsea_neg_score",
                           "gsea_score",
                           "gsea_pos_pval",
                           "gsea_neg_pval",
                           "gsea_combined_pval",
                           "padj")

        # FILTER AND SORT
        filtered_result <- result_padj[, relevant_cols] %>% arrange(-gsea_score)

        new_cols <- c("AS Signature",
                      "Positive ES",
                      "Negative ES",
                      "SPARKS ES",
                      "Positive p-val",
                      "Negative p-val",
                      "Combined p-val",
                      "FWER")
        colnames(filtered_result) <- new_cols

        return(filtered_result)
    })

    # process manual colors
    select_gene_colors <- reactive({
        color_vector <- rep("red", length(input$custom_highlight_rbps))
        names(color_vector) <- input$custom_highlight_rbps

        return(color_vector)
    })

    output$barplot <- renderPlot({
        req(input$num_plot)
        generate_enrichment_barplot(input_sparks_result(),
                                    num_plot = input$num_plot,
                                    select_genes = input$custom_highlight_rbps,
                                    manual_colors = select_gene_colors(),
                                    text_scale_factor = 1.5,
                                    select_gene_marker = T)
    })


}

# Run the application
shinyApp(ui = ui, server = server)

