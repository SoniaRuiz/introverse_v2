#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)


source("backend.R")


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  tags$style(HTML('.shiny-split-layout>div { overflow:visible; }')), 
  
   # Application title
   titlePanel("IntroVerse v2"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        p(strong("Insert the coordinates of an annotated intron, alternative 5' or alternative 3' splicing event of interest:")),
        shiny::splitLayout(id="chr_strand_tab1",
            shiny::selectInput(inputId = "chr_tab1",
                               label = "Chr",
                               selected = 1,
                               choices = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"),
                               multiple = F),
            shiny::selectInput(inputId = "strand_tab1",
                               label = "Strand",
                               choices = c("+", "-"),
                               selected = "-",
                               multiple = F)
        ),
        shiny::splitLayout(id="start_end_tab1",
           numericInput(inputId = "start_tab1",
                        label = "Start",
                        value = 155235308),
           numericInput(inputId = "end_tab1",
                        label = "End",
                        value = 155235680)
        ),
        hr(),
        p(strong("Select the database to retrieve splicing information about the event indicated above:")),
        
        shiny::splitLayout(id="database_selection_all_tab1",
                           shiny::radioButtons(inputId = "database_select",
                                               choices = c("All",
                                                           "GTEx v8",
                                                           "GTEx v8 Age",
                                                           "ENCODE shRNA",
                                                           "Alzheimer's Disease/Control"),
                                               label = "Database",
                                               selected = "ENCODE shRNA")
                           ),
    
        actionButton(inputId = "acceptButton_tab1", label = "Accept"),
        
        ## Hidden inputs for the transcript visualisation
        shinyjs::hidden(
          p(id = "junID_tab1", ""),
          p(id = "transcriptENS_tab1", ""),
          p(id = "databaseName_tab1", "")
        )
        
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         uiOutput("queryResults") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
         
         shinyBS::bsModal(id = "modalVisualiseTranscript_tab1", 
                          title = NULL,
                          trigger = NULL,
                          size = "large",
                          plotOutput("modalVisualiseTranscript_tab1") %>% 
                            shinycssloaders::withSpinner(color="#0dc5c1"))
                          #downloadButton(outputId = 'downloadPlot', label = 'Download (.svg)'))
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  #updateSelectizeInput(session, 'gtex_table_list', choices = gtex_table_list, server = TRUE, selected = gtex_table_list[1])
  #updateSelectizeInput(session, 'encode_table_list', choices = encode_table_list, server = TRUE, selected = encode_table_list[1])
  
  
  ##############################################################################
  ## MODAL POP-UPS
  ##############################################################################
  
  visualiseTranscriptPlot <- function() {
    visualise_transcript(junID = str_replace_all(string = input$junID_tab1, pattern = "%20", replacement = " "),
                         transcript_ENS = str_replace_all(string = input$transcriptENS_tab1, pattern = "%20", replacement = " "),
                         database_name = str_replace_all(string = input$databaseName_tab1, pattern = "%20", replacement = " "))
  }
  
  output$modalVisualiseTranscript_tab1 <- renderPlot({
    shinyjs::runjs(code = '$("#modalVisualiseTranscript_tab1 > img").addClass("d-none");')
    shinyjs::runjs(code = '$("#modalVisualiseTranscript_tab1 > img").remove();')
    visualiseTranscriptPlot()
  }, width = "auto", height = "auto")
  
  
  ##############################################################################
  ## RETRIEVE INFO FROM DATABASES
  ##############################################################################
  
  toListen <- eventReactive(list(input$acceptButton_tab1), {
    

    
    query_results <- query_database(input$chr_tab1, 
                                    input$start_tab1, 
                                    input$end_tab1, 
                                    input$strand_tab1,
                                    input$database_select)
    
    
    if ( query_results %>% nrow > 0 ) {
      
      
      ## Info about the junction that is common to all databases
      v <- setup_ui_details_section(query_results)
      
      i <- (v %>% length()) + 1
      
      v[[i]] <- div(
        br(),
        hr(),
        h3("Results"),
        br()
      )
      i <- i + 1
      
      ## Info particular to each database
      v <- setup_ui_results_section(query_results, v, i)
      
     # print(v)
      v
      
    } else {
      p(strong("No results found!"))
    }
    
    
    
  })
    
  
  output$queryResults <- renderUI({
    toListen()
  })

  
}

# Run the application 
shinyApp(ui = ui, server = server)

