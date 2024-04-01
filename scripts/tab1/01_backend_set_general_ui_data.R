
set_UI_distance_information <- function(jxn_info) {
  
  if ( str_detect(string = jxn_info$junID_type %>% unique(), pattern = "novel" ) ) {
    
    info_jxn_header <- c("Distance to annotated intron: ", "Reference intron: ")
    info_jxn <- c(paste0(jxn_info$distance %>% unique, " bp."),
                  paste0(paste(jxn_info$ref_coordinates %>% unique, collapse = " ")))
    if (jxn_info$junID_type %>% unique() == "novel_donor") {
      text_color = "#35B779FF"
    } else if ( jxn_info$junID_type %>% unique() == "novel_acceptor") {
      text_color = "#64037d"
    }
    
  } else {
    info_jxn_header  = c("Mis-spliced?","")
    if ( any(jxn_info$ref_type %>% unique == "never" )) {
      info_jxn <- c("No (accurately spliced intron)","")
      text_color = c("blue", "green")
    } else {
      if ( any(jxn_info$ref_type %>% unique == "both" )) {
        info_jxn <- c(paste0("Yes (both splice sites)."),"") 
      } else {
        info_jxn <- c(paste0("Yes (",jxn_info$ref_type %>% unique, " splice sites)."),"")  
      }
      text_color = c("blue", "orange")
    }
    
  } 
  
  return(data.frame("distance_label" = info_jxn_header,
                    "distance_info" = info_jxn,
                    "color" = text_color))
  
}

set_UI_transcript_information <- function(query_results) {
  
  
  ## Get transcript info
  transcript_span_tag <- map(query_results$transcript_ENS %>% unique %>% sort(), function(transcript) {
    ## TODO if transcript has exons, utr and cds
    span(tags$a(href = "#", 
                onclick = "var ID_jxn= $('#span_junID')[0].innerText; 
                  Shiny.setInputValue('junID_tab1', ID_jxn);
                  var ID_database= $('#h4_databaseID')[0].innerText; 
                  Shiny.setInputValue('databaseName_tab1', ID_database);
                  Shiny.setInputValue('transcriptENS_tab1', this.text );
                  var junType= $('#span_junType')[0].innerText; 
                  Shiny.setInputValue('junType_tab1', junType);
                  $('#modalVisualiseTranscript_tab1').modal('show');",
                transcript ), HTML('&nbsp;'))
    
    
  })
  return(transcript_span_tag)
  
}

set_UI_clip_information <- function(query_results) {
  
  ## Get binding sites info
  encori_overlaps <- NULL
  encori_span_tag <- NULL

  ## Connect to the CLIP data database and retrieve data from the current gene
  database_path <- file.path("./database/clip_data.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  if ( DBI::dbExistsTable(conn = con, name = query_results$gene_name %>% unique) ) {
    
    query <- paste0("SELECT * FROM '", query_results$gene_name %>% unique, "'")
    db_clip_data_gr <- DBI::dbGetQuery(con, query) %>%
      GenomicRanges::GRanges()
    DBI::dbDisconnect(conn = con)
    
    if ( db_clip_data_gr %>% length() > 0 ) {
      encori_overlaps <- GenomicRanges::findOverlaps(query = query_results %>% 
                                                       drop_na(ref_junID) %>% 
                                                       distinct(coordinates, .keep_all = T) %>% 
                                                       dplyr::select(seqnames, start, end, strand ) %>%
                                                       GenomicRanges::GRanges(),
                                                     subject = db_clip_data_gr,
                                                     maxgap = 100,
                                                     type = "any")
      
      if ( encori_overlaps %>% length() > 0 ) {
        
        encori_hits <- db_clip_data_gr[S4Vectors::subjectHits(encori_overlaps),]$RBP %>% unique
        
        encori_span_tag <- div(p(strong("Contains binding sites for:"), HTML('&nbsp;'),
                                 span(tags$a(href = "#", 
                                             onclick = "
                                       var ID_jxn= $('#span_junID')[0].innerText; 
                                       Shiny.setInputValue('junID_tab1', ID_jxn);
                                       var gene_Name= $('#span_geneName')[0].innerText; 
                                       Shiny.setInputValue('geneName_tab1', gene_Name);
                                       var transcript_Name= $('#transcript_list span a')[0].innerText; 
                                       Shiny.setInputValue('transcriptENS_tab1', transcript_Name );
                                       var junType= $('#span_junType')[0].innerText; 
                                       Shiny.setInputValue('junType_tab1', junType);
                                       $('#modalVisualiseCLIP_tab1').css({'height':'90vh'});
                                       $('#modalVisualiseCLIP_tab1').modal('show');",
                                             paste0(encori_hits %>% length(), " RBPs") ), HTML('&nbsp;'))))
      } 
    } 
  } 
  
  if ( is.null(encori_span_tag) ) {
    encori_span_tag <- div(p(strong("Contains binding sites for:"), HTML('&nbsp;'), "No CLIP data available for this gene."))
  }
  
  return(encori_span_tag)
}

set_UI_clinvar_information <- function(jxn_info) {
  
  
  ## Get clinvar info
  clinvar_span_tag <- NULL
  
  if (jxn_info %>% drop_na(clinvar) %>% pull(clinvar) %>% unique == 1) {
    
    #print(jxn_info$clinvar)
    #print(jxn_info$clinvar_locus)
    
    if ( any(jxn_info %>% names == "clinvar_locus") ) {
      
      clinvar_hits <- jxn_info %>% drop_na("clinvar_locus") %>% pull(clinvar_locus) %>% unique
      
      clinvar_span_tag <-div(p(strong("ClinVar variants:"), HTML('&nbsp;'),
                               span(tags$a(href = "#", 
                                           onclick = "var ID_jxn= $('#span_junID')[0].innerText; 
                                           Shiny.setInputValue('junID_tab1', ID_jxn);
                                           var ID_database= $('#h4_databaseID')[0].innerText; 
                                           Shiny.setInputValue('databaseName_tab1', ID_database);
                                           Shiny.setInputValue('clinvarlocus_tab1', this.text );
                                           var junType= $('#span_junType')[0].innerText; 
                                           Shiny.setInputValue('junType_tab1', junType);
                                           var gene_Name= $('#span_geneName')[0].innerText; 
                                           Shiny.setInputValue('geneName_tab1', gene_Name);
                                           $('#modalVisualiseClinVar_tab1').modal('show');",
                                           paste0(clinvar_hits)), HTML('&nbsp;'))))
    }
  }
  
  if ( is.null(clinvar_span_tag) ) {
    clinvar_span_tag <- div(p(strong("ClinVar variants:"), HTML('&nbsp;'), "No ClinVar data."))
  }
  
  return(clinvar_span_tag)
}