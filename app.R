#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)

source("scripts/tab1/99_backend_utils.R")
source("backend.R")
source("scripts/tab1/01_backend_set_general_ui_data.R")
source("scripts/tab1/02_backend_visualise_in_transcript.R")
source("scripts/tab1/03_backend_download_data.R")

source("scripts/tab2/01_backend_tab2.R")



# Define UI for application that draws a histogram
ui <- navbarPage("IntroVerse V2",
                 inverse=TRUE,
                 id = "introverse_menu",
                 useShinyjs(),
                 
                 tags$head(
                   tags$link(rel = "stylesheet", type = "text/css", href = "api.css")
                 ),
                 
                 
                 tabPanel("Basic Search",
                          
                          
                          shinybusy:: use_busy_spinner(
                            spin_id = "main_spinner",
                            spin = "circle",
                            color = "#0dc5c1",
                            position = "full-page"
                          ),
                          
                          # Sidebar with a slider input for number of bins 
                          sidebarLayout(
                            sidebarPanel(
                              width = 4,
                              p(strong("Insert the coordinates of an annotated intron, alternative 5' or alternative 3' splicing event of interest (hg38):")),
                              shiny::splitLayout(id="chr_strand_tab1",
                                                 # shiny::textInput(inputId = "jxn_coordinates_tab1",
                                                 #                  label = "Coordinates (hg38):",
                                                 #                  value = "chr16:86478133-86508654:-")
                                                 shiny::selectInput(inputId = "chr_tab1",
                                                                    label = "Chr",
                                                                    selected = 1,
                                                                    choices = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"),
                                                                    multiple = F),
                                                 numericInput(inputId = "start_tab1",
                                                              label = "Start",
                                                              value = 154597977),
                                                 numericInput(inputId = "end_tab1",
                                                              label = "End",
                                                              value = 154598401),
                                                 shiny::selectInput(inputId = "strand_tab1",
                                                                    label = "Strand",
                                                                    choices = c("+", "-"),
                                                                    selected = "-",
                                                                    multiple = F)
                              #),
                              #shiny::splitLayout(id="start_end_tab1",
                                                 
                              ),
                              hr(),
                              p(strong("Select the database to retrieve splicing information about the splicing event indicated above.")),
                              
                              shiny::splitLayout(id="database_selection_all_tab1",
                                                 shiny::radioButtons(inputId = "database_select",
                                                                     choices = c("All",
                                                                                 "TCGA",
                                                                                 "GTEx v8",
                                                                                 "GTEx v8 Age",
                                                                                 "ENCODE shRNA",
                                                                                 "Alzheimer's Disease/Control",
                                                                                 "Parkinson's Disease/Control",
                                                                                 "LRRK2 G2019S mutation"),
                                                                     label = "Database:",
                                                                     selected = "TCGA")
                              ),
                              hr(),
                              actionButton(inputId = "acceptButton_tab1", class = "btn-primary", label = "Accept"),
                              
                              ## Hidden inputs for the transcript visualisation
                              shinyjs::hidden(
                                p(id = "junID_tab1", ""),
                                p(id = "geneName_tab1", ""),
                                p(id = "junType_tab1", ""),
                                p(id = "transcriptENS_tab1", ""),
                                p(id = "databaseName_tab1", ""),
                                p(id = "clinvarlocus_tab1", ""),
                                p(id = "phastCons17_5ss_tab1",""),
                                p(id = "phastCons17_3ss_tab1",""),
                                p(id = "CDTS_5ss_tab1",""),
                                p(id = "CDTS_3ss_tab1",""),
                                p(id = "MES_5ss_tab1",""),
                                p(id = "MES_3ss_tab1",""),
                                p(id = "donor_sequence_tab1",""),
                                p(id = "acceptor_sequence_tab1","")
                              )
                              
                            ),
                            
                            # Show the output
                            mainPanel(
                              width = 8,
                              shiny::uiOutput(outputId = "queryResults"), 
                              tabsetPanel(
                                id = "switcher",
                                type = "tabs" 
                              ),
                              
                              
                              shinyBS::bsModal(id = "modalVisualiseTranscript_tab1", 
                                               title = "Visualise the splicing event within the selected transcript",
                                               trigger = NULL,
                                               size = "large",
                                               
                                               shinybusy:: use_busy_spinner(
                                                 spin_id = "modal_spinner",
                                                 spin = "circle",
                                                 color = "#0dc5c1",
                                                 position = "full-page"
                                               ),
                                               
                                               plotOutput("modalVisualiseTranscript_tab1", width = "auto", height = "70vh")),
                              
                              shinyBS::bsModal(id = "modalVisualiseCLIP_tab1", 
                                               title = "RBP-RNA interactions supported by CLIP-seq data",
                                               
                                               trigger = NULL,
                                               size = "large",
                                               
                                               shinybusy:: use_busy_spinner(
                                                 spin_id = "modal_spinner",
                                                 spin = "circle",
                                                 color = "#0dc5c1",
                                                 position = "full-page"
                                               ),
                                               
                                               
                                               plotOutput("modalVisualiseCLIP_tab1", width = "auto", height = "70vh"),
                                               footer = downloadButton(outputId = 'downloadCLIPData', label = 'Download all regulatory data',
                                                                       class = "btn-primary")
                                               ),
                              
                              shinyBS::bsModal(id = "modalVisualiseClinVar_tab1", 
                                               title = "ClinVar pathogenic/likely-pathogenic splicing variants",
                                               
                                               trigger = NULL,
                                               size = "large",
                                               shinybusy:: use_busy_spinner(
                                                 spin_id = "modal_spinner",
                                                 spin = "circle",
                                                 color = "#0dc5c1",
                                                 position = "full-page"
                                               ),
                                               
                                               plotOutput("modalVisualiseClinVar_tab1", width = "auto", height = "65vh")),
                              
                              shinyBS::bsModal(id = "modalVisualiseMORE_tab1", 
                                               title = "More information about the splicing event selected",
                                               trigger = NULL,
                                               size = "large",
                                               uiOutput(outputId = "modalVisualiseMORE_tab1", width = "auto", height = "70vh") )
                             
                            )
                          )
                          ),
                 tabPanel("Datasets",
                          
                          shinybusy:: use_busy_spinner(
                            spin_id = "main_spinner_tab2",
                            spin = "circle",
                            color = "#0dc5c1",
                            position = "full-page"
                          ),
                          
                          fluidPage(
                            navlistPanel(
                              id = "tab_panel_DB_info",
                              "Databases",
                              tabPanel(
                                title = "TCGA",
                                fluidRow(
                                  column(9,
                                         h1("TCGA"),
                                         h3(strong("The Cancer Genome Atlas Program")),
                                         p("The Cancer Genome Atlas (TCGA) is a landmark cancer genomics program, molecularly characterized over 20,000 primary cancer and matched normal samples spanning 33 cancer types"),
                                         plotOutput("tcga_metadata_output_tab2", width = "auto", height = "50vh"),
                                         
                                         uiOutput(outputId = "tcga_cancer_types_output_tab2")
                                  ))),
                              tabPanel(
                                title = "GTEx v8",
                                fluidRow(
                                  column(9,
                                         h1("GTEx v8"),
                                         h3(strong("Genotype Tissue Expression Project v8")),
                                         p("GTEx project collected samples from up to 54 non-diseased tissue sites across nearly 1,000 deceased individuals. Gene expression of each tissue was assessed by RNA sequencing (bulk RNA-seq)."),
                                         plotOutput("gtex_metadata_output_tab2", width = "auto", height = "50vh")
                                  ))),
                              tabPanel(
                                title = "GTEx v8 Age",
                                fluidRow(
                                  column(9,
                                         h1("GTEx v8 Age"),
                                         h3(strong("Genotype Tissue Expression Project v8. ")),
                                         p("GTEx samples have been clustered by age in three age supergroups: '20-39', '40-59' and '60-79' years-old. "),
                                         plotOutput("gtex_age_metadata_output_tab2", width = "auto", height = "50vh")
                                  ))),
                              tabPanel(
                                title = "ENCODE shRNA",
                                fluidRow(
                                  column(9,
                                         h1("ENCODE shRNA"),
                                         h3(strong("shRNA knockdown followed by RNA-seq (shRNA RNA-seq) experiments")),
                                         p("RNA-seq on K562 and HepG2 cells treated with an shRNA knockdown against different RNA-Binding Proteins (RBPs) as well as against no target and followed by RNA-seq (ENCODE platform)."),
                                         plotOutput("ENCODE_metadata_output_tab2", width = "auto", height = "50vh")
                                  ))),
                              tabPanel(
                                title = "Alzheimer's Disease/Control",
                                fluidRow(
                                  column(9,
                                         h1("Alzheimer's Disease/Control"),
                                         h3(strong("NCBI: "),"SRP10094"),
                                         p("Bulk RNA-sequencing was purified from fusiform gyrus tissue sections of autopsy-confirmed Alzheimer's cases and neurologically normal age-matched controls."),
                                         plotOutput("SRP100948_metadata_output_tab2", width = "auto", height = "50vh")
                                  ))),
                              tabPanel(
                                title = "Parkinson's Disease/Control",
                                fluidRow(
                                  column(9,
                                         h1("Parkinson's Disease/Control"),
                                         h3(strong("NCBI: "),"SRP058181"),
                                         p("Bulk RNA-sequencing expression of human post-mortem BA9 brain tissue for Parkinson Disease and neurologically normal individuals."),
                                         plotOutput("SRP058181_metadata_output_tab2", width = "auto", height = "50vh")
                                  ))),
                              tabPanel(
                                title = "LRRK2 G2019S mutation",
                                fluidRow(
                                  column(9,
                                         h1("LRRK2 G2019S mutation"),
                                         h3(strong("NCBI:"),"SRP151040"),
                                         #p(strong("Abstract")),
                                         p("Short-read RNA-sequencing of induced pluripotent stem cells (iPSC) and iPSC-derived astrocytes from control and Parkinson's disease patients carrying LRRK2 G2019S point mutation"),
                                         #p(strong("Overall study design"), "Human iPSC-derived astrocyte identity was obtained by comparing human iPSC from the same control and Parkinson's disease LRRK2 G2019S patients and human commercial astrocytes. mRNA profiles were generated by deep sequencing using Illumina HiSeq2000"),
                                         #p(strong("Original study goal: "), 
                                        #   "The goal of this study is to compare the NGS-derived from transcriptome profiling (RNA-seq) of human iPSC, human iPSC-derived astrocytes from control and 
                                        #   Parkinson's disease LRRK2 G2019S, and human commercial astrocytes to gain insight into the identity of human iPSC-derived astrocytes in vitro during the 
                                        #   differentiation process."),
                                         #p("Total RNA was assayed for quantity and quality using Qubit RNA HS Assay (Life Technologies) and RNA 6000 Nano Assay on a Bioanalyzer 2100. 
                                         #The RNASeq libraries were prepared from total RNA using the TruSeq®Stranded mRNA LT Sample Prep Kit (Illumina Inc., Rev.E, October 2013).  
                                         #Libraries were sequenced on HiSeq2000 (Illumina, Inc) in paired-end mode with a read length of 2x76bp using TruSeq SBS Kit v4. 
                                         #Over 30 million paired-end reads for each sample in a fraction of a sequencing v4 flow cell lane, following the manufacturer's protocol."), 
                                         #p(strong("Study Conclusion:"), "These results suggest that iPSC-derived astrocytes resemble human commercial astrocytes validating the differentiating protocol used."),
                                           
                                         plotOutput("SRP151040_metadata_output_tab2", width = "auto", height = "50vh")
                                  )))
                              
                            ))
                 )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  #updateSelectizeInput(session, 'gtex_table_list', choices = gtex_table_list, server = TRUE, selected = gtex_table_list[1])
  #updateSelectizeInput(session, 'encode_table_list', choices = encode_table_list, server = TRUE, selected = encode_table_list[1])
  
  
  shiny::updateTabsetPanel(inputId = "introverse_menu", selected = "Basic Search")
  
  ##############################################################################
  ## MODAL POP-UP 
  ##############################################################################
  
  ## TRANSCRIPT VISUALISATION ----------------------------------------------------
  
  visualiseTranscriptPlot <- function() {
    transcript_plot <- visualise_transcript(junID = str_replace_all(string = input$junID_tab1, pattern = "%20", replacement = " "),
                                            transcript_ENS = str_replace_all(string = input$transcriptENS_tab1, pattern = "%20", replacement = " "),
                                            database_name = str_replace_all(string = input$databaseName_tab1, pattern = "%20", replacement = " "),
                                            junType = str_replace_all(string = input$junType_tab1, pattern = "%20", replacement = " "))
    shinybusy::hide_spinner(spin_id = "modal_spinner")
    return(transcript_plot)
  }
  
  output$modalVisualiseTranscript_tab1 <- renderPlot({
    shinybusy::show_spinner(spin_id = "modal_spinner")
    shinyjs::runjs(code = '$("#modalVisualiseTranscript_tab1 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#modalVisualiseTranscript_tab1 > img").remove();')
    visualiseTranscriptPlot()
  }, res = 75)
  
  
  
  ## CLIP VISUALISATION ---------------------------------------------------------
  
  visualiseCLIP_Plot <- function() {
    CLIP_plot <- visualise_CLIP(junID = str_replace_all(string = input$junID_tab1, pattern = "%20", replacement = " "),
                                transcript_ENS = str_replace_all(string = input$transcriptENS_tab1, pattern = "%20", replacement = " "),
                                geneName = str_replace_all(string = input$geneName_tab1, pattern = "%20", replacement = " "),
                                junType = str_replace_all(string = input$junType_tab1, pattern = "%20", replacement = " "))
    shinybusy::hide_spinner(spin_id = "modal_spinner")
    return(CLIP_plot)
  }
  
  output$modalVisualiseCLIP_tab1 <- renderPlot({
    shinybusy::show_spinner(spin_id = "modal_spinner")
    shinyjs::runjs(code = '$("#modalVisualiseCLIP_tab1 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#modalVisualiseCLIP_tab1 > img").remove();')
    visualiseCLIP_Plot()
  }, res = 75)
  
  
  ## CLINVAR VISUALISATION ----------------------------------------------------
  
  visualiseClinVarPlot <- function() {
    clinvar_plot <- visualise_clinvar(junID = str_replace_all(string = input$junID_tab1, pattern = "%20", replacement = " "),
                                      clinvar_locus = str_replace_all(string = input$clinvarlocus_tab1, pattern = "%20", replacement = " "),
                                      geneName = str_replace_all(string = input$geneName_tab1, pattern = "%20", replacement = " "),
                                      database_name = str_replace_all(string = input$databaseName_tab1, pattern = "%20", replacement = " "),
                                      junType = str_replace_all(string = input$junType_tab1, pattern = "%20", replacement = " "))
    shinybusy::hide_spinner(spin_id = "modal_spinner")
    return(clinvar_plot)
  }
  
  output$modalVisualiseClinVar_tab1 <- renderPlot({
    shinybusy::show_spinner(spin_id = "modal_spinner")
    shinyjs::runjs(code = '$("#modalVisualiseClinVar_tab1 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#modalVisualiseClinVar_tab1 > img").remove();')
    visualiseClinVarPlot()
  }, res = 72)
  
  
  ## 'MORE...' VISUALISATION ----------------------------------------------------
  
  visualiseMORE_Table <- function() {
    
    DT::datatable(data = data.frame("Property" = c("PhastCons17 5'ss (100 bp)","PhastCons17 3'ss (100 bp)","CDTS 5'ss (100 bp)",
                                                   "CDTS 3'ss (100 bp)","MaxEntScan 5'ss (9 bp)","MaxEntScan 3'ss (23 bp)", 
                                                   "5'ss sequence (3bp | 6bp)", "3'ss sequence (20bp | 3bp)" ),
                                    "Value" = c(input$phastCons17_5ss_tab1,input$phastCons17_3ss_tab1,input$CDTS_5ss_tab1,
                                                input$CDTS_3ss_tab1,input$MES_5ss_tab1,input$MES_3ss_tab1,
                                                input$donor_sequence_tab1,input$acceptor_sequence_tab1),
                                    "Description" = c("Inter-species conservation across 17 primates. Calculated using the 100 bp sequence overlapping the 5'ss junction to 100 bp downstream the intron.\n
                                                      The higher is the PhastCons17 score, the more conserved is the 100 bp sequence evaluated across primates.",
                                                      
                                                      "Inter-species conservation across 17 primates. Calculated using the 100 bp sequence overlapping the 3'ss junction to 100 bp upstream the intron.",
                                                      
                                                      "Sequence constraint across humans. Calculated using the 100 bp sequence overlapping the 5'ss junction to 100 bp downstream the intron.\n
                                                      The higher is the CDTS score, the more sequence variation has been found across humans for the 100 bp sequence evaluated.",
                                                      
                                                      "Sequence constraint across humans. Calculated using the 100 bp sequence overlapping the 3'ss junction to 100 bp upstream the intron.",
                                                      
                                                      "Sequence motif similarity between the 9 bp sequence overlapping the 5'ss splice sequence and annotated splice sites. 
                                                      The higher is the score, the higher is the likelihood of the 5' splice site being recognised by the spliceosome machinery.",
                                                      
                                                      "Sequence motif similarity between the 23 bp sequence overlapping the 3'ss splice sequence and annotated splice sites.", 
                                                      
                                                      "9 bp motif DNA sequences located at the 5'ss (-3/+6 bp). Obtained using the Human Primary DNA Assembly hg38.",
                                                      "23 bp motif DNA sequences located at the 3'ss (-20/+3 bp). Obtained using the Human Primary DNA Assembly hg38." )),
                  rownames = F, 
                  extensions = c('Buttons', 'Responsive'),
                  
                  options = list(
                    pageLength = 10,
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel')
                  )) %>%
      return()
  }
  
  output$modalVisualiseMORE_tab1 <- renderUI({
    shinyjs::runjs(code = '$("#modalVisualiseMORE_tab1 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#modalVisualiseMORE_tab1 > img").remove();')
    visualiseMORE_Table() %>% return()
  })

  
  ##############################################################################
  ## RETRIEVE SPLICING INFO FROM DATABASES
  ##############################################################################

  toListen <- eventReactive(list(input$acceptButton_tab1), {
    
    ## Query databases
    query_results <- query_database(chr = input$chr_tab1,
                                    start = input$start_tab1,
                                    end = input$end_tab1,
                                    strand = input$strand_tab1,
                                    input$database_select)
    
    if ( query_results %>% nrow > 0 ) {
    
      ## Info about the junction that is common to all databases
      v <- setup_UI_details_section(query_results)
      i <- (v %>% length()) + 1
      v[[i]] <- div(
        br(),
        hr(),
        h3("Splicing information per database"),
        br()
      )
      
      ## Create a tab per database containg the splicing info particular to each database
      result_tabs <- setup_UI_results_section(query_results)
      
      if ( result_tabs %>% length > 0 ) {
        
        lapply(1:length(result_tabs), FUN =  function(j) {
          shiny::insertTab(
            inputId = "switcher",
            result_tabs[[j]],
            target = NULL,
            select = T
          )
        })
        
        ## Show new tab results
        shinyjs::show(id = "switcher")
        shinyjs::show(selector = ".tabbable div.tab-content")
      } 
      
      ## Hide loading spinner
      shinybusy::hide_spinner(spin_id = "main_spinner")
      
      ## Return results
      return(v)
      
    } else {
      shinybusy::hide_spinner(spin_id = "main_spinner")
      p(strong("No results found!"))
    }
    
  })

  output$queryResults <- renderUI({
    
    ## Remove previous tabs
    for (tab_name in database_equivalences$name){
      shiny::removeTab(inputId = "switcher", target = tab_name)  
    }
    
    
    ## Hide previous results
    shinyjs::hide(id = "switcher")
    shinyjs::hide(selector = ".tabbable div.tab-content")
    
    ## Show loading spinner
    shinybusy::show_spinner(spin_id = "main_spinner")
    
    ## Query database
    toListen()
    
    
  })
  
  ###########################################
  ## RETRIEVE METADATA OUTPUT
  ###########################################
  
  toListen_getDBMetadata <- eventReactive(list(input$tab_panel_DB_info), {
    if (input$tab_panel_DB_info == "GTEx v8") {
      data_to_plot <- visualise_metadata_sample_clusters_gtex(input$tab_panel_DB_info)    
    } else {
      data_to_plot <- visualise_metadata_sample_clusters(input$tab_panel_DB_info)  
    }
    
    shinybusy::hide_spinner(spin_id = "main_spinner_tab2")
    data_to_plot %>% return()
  })
  
  output$tcga_metadata_output_tab2 <- renderPlot({
    shinybusy::show_spinner(spin_id = "main_spinner_tab2")
    shinyjs::runjs(code = '$("#tcga_metadata_output_tab2 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#tcga_metadata_output_tab2 > img").remove();')
    toListen_getDBMetadata()
  })
  output$tcga_cancer_types_output_tab2 <- renderUI({
    div(h3(strong("Cancer types")),
        DT::datatable(data = table_cancer_types(),
                      rownames = F)) %>% return()
  })
  
  
  output$gtex_metadata_output_tab2 <- renderPlot({
    shinybusy::show_spinner(spin_id = "main_spinner_tab2")
    shinyjs::runjs(code = '$("#gtex_metadata_output_tab2 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#gtex_metadata_output_tab2 > img").remove();')
    toListen_getDBMetadata()
  })
  
  output$gtex_age_metadata_output_tab2 <- renderPlot({
    shinybusy::show_spinner(spin_id = "main_spinner_tab2")
    shinyjs::runjs(code = '$("#gtex_age_metadata_output_tab2 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#gtex_age_metadata_output_tab2 > img").remove();')
    toListen_getDBMetadata()
  })
  
  output$ENCODE_metadata_output_tab2 <- renderPlot({
    shinybusy::show_spinner(spin_id = "main_spinner_tab2")
    shinyjs::runjs(code = '$("#ENCODE_metadata_output_tab2 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#ENCODE_metadata_output_tab2 > img").remove();')
    toListen_getDBMetadata()
  })
  
  output$SRP100948_metadata_output_tab2 <- renderPlot({
    shinybusy::show_spinner(spin_id = "main_spinner_tab2")
    shinyjs::runjs(code = '$("#SRP100948_metadata_output_tab2 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#SRP100948_metadata_output_tab2 > img").remove();')
    toListen_getDBMetadata()
  })
  
  output$SRP058181_metadata_output_tab2 <- renderPlot({
    shinybusy::show_spinner(spin_id = "main_spinner_tab2")
    shinyjs::runjs(code = '$("#SRP058181_metadata_output_tab2 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#SRP058181_metadata_output_tab2 > img").remove();')
    toListen_getDBMetadata()
  })
  
  output$SRP151040_metadata_output_tab2 <- renderPlot({
    shinybusy::show_spinner(spin_id = "main_spinner_tab2")
    shinyjs::runjs(code = '$("#SRP151040_metadata_output_tab2 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#SRP151040_metadata_output_tab2 > img").remove();')
    toListen_getDBMetadata()
  })
  
  
  ###########################################
  ## DOWNLOAD HANDLERS
  ###########################################

  output$downloadCLIPData <- downloadHandler(
    filename = function() {
      paste("CLIPdata_",input$junID_tab1, "_", input$geneName_tab1,".csv", collapse = "")
    },
    content = function(file) {
      CLIP_data <- download_CLIP_data(input$junID_tab1, input$geneName_tab1)
      write.csv(x = CLIP_data, 
                file, 
                row.names = F,
                fileEncoding = "UTF-8")
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

