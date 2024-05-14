


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
    if ( all(jxn_info$ref_type %>% unique == "never" ) ) {
      info_jxn <- c("No (accurately spliced intron)","")
      text_color = c("blue", "green")
    } else {
      if ( any(jxn_info$ref_type %>% unique == "both" ) ) {
        info_jxn <- c(paste0("Yes (both splice sites)."),"") 
      } else {
        splice_sites <- jxn_info$ref_type %>% unique
        if (any(splice_sites == "never")) {
          splice_sites[which(splice_sites=="never")] <- "accurate"
        }
        info_jxn <- c(paste0("Yes (", paste(splice_sites, collapse = " and "), " splicing)."),"")  
      }
      text_color = c("blue", "orange")
    }
    
  } 
  
  return(data.frame("distance_label" = info_jxn_header,
                    "distance_info" = info_jxn,
                    "color" = text_color))
  
}

set_UI_transcript_information <- function(query_results) {
  
  query_results_transcript <- query_results %>%
    dplyr::select(transcript_ENS, MANE) %>% 
    distinct(transcript_ENS, .keep_all = T) %>% 
    arrange(desc(MANE), transcript_ENS) 
  
  transcript_id <- ""
  MANE_text = ""
  
  if ( any(query_results_transcript$MANE == 1) ) {
    transcript_id <- query_results_transcript %>%
      filter(MANE == 1) %>%
      pull(transcript_ENS)
    MANE_text = " (MANE)"
  } else {
    transcript_id <- (query_results_transcript %>%
      pull(transcript_ENS))[1]
  }
  
  query_results_transcript <- query_results_transcript %>%
    filter(transcript_ENS != transcript_id)
  
  
  transcript_span_tag <- 
    span(
      tags$a(id = "atag_transcriptID",
             href = "#", 
      onclick = "var ID_jxn= $('#span_junID')[0].innerText; 
      Shiny.setInputValue('junID_tab1', ID_jxn);
      var ID_database= $('#h4_databaseID')[0].innerText; 
      Shiny.setInputValue('databaseName_tab1', ID_database);
      Shiny.setInputValue('transcriptENS_tab1', this.text );
      var junType= $('#span_junType')[0].innerText; 
      Shiny.setInputValue('junType_tab1', junType);
      Shiny.setInputValue('transcriptMANE_tab1', 1);
      $('#modalVisualiseTranscript_tab1').modal('show');", transcript_id), 
      MANE_text)
    
  
  transcript_span_tag <- span(transcript_span_tag, "and",
                             span(tags$a(href = "#", 
                             onclick = "var ID_jxn= $('#span_junID')[0].innerText; 
                             Shiny.setInputValue('junID_tab1', ID_jxn);
                             var ID_database= $('#h4_databaseID')[0].innerText; 
                             Shiny.setInputValue('databaseName_tab1', ID_database);
                             var transcript_ID = $('#atag_transcriptID')[0].innerText;
                             Shiny.setInputValue('transcriptENS_tab1', transcript_ID );
                             var junType= $('#span_junType')[0].innerText; 
                             Shiny.setInputValue('junType_tab1', junType);
                             var gene_Name= $('#span_geneName')[0].innerText; 
                             Shiny.setInputValue('geneName_tab1', gene_Name);
                             Shiny.setInputValue('transcriptMANE_tab1', 0);
                             $('#modalVisualiseTranscript_tab1').modal('show');", query_results_transcript %>% filter(MANE == 0) %>% nrow(), " more transcripts"#, HTML('&nbsp;')
                             )))
  
  return(transcript_span_tag)
  
  
  ## Get transcript info
  # transcript_span_tag <- map(.x = split(query_results_transcript, 1:nrow(query_results_transcript)), 
  #                            .f = function(row) {
  #   
  #   transcript_id <- row$transcript_ENS
  #   
  #   MANE_text <- ""
  #   if (row$MANE) { MANE_text <- " (MANE)" }
  #   
  #   span(tags$a(href = "#", 
  #               onclick = "var ID_jxn= $('#span_junID')[0].innerText; 
  #                 Shiny.setInputValue('junID_tab1', ID_jxn);
  #                 var ID_database= $('#h4_databaseID')[0].innerText; 
  #                 Shiny.setInputValue('databaseName_tab1', ID_database);
  #                 Shiny.setInputValue('transcriptENS_tab1', this.text );
  #                 var junType= $('#span_junType')[0].innerText; 
  #                 Shiny.setInputValue('junType_tab1', junType);
  #                 $('#modalVisualiseTranscript_tab1').modal('show');",
  #               transcript_id ), MANE_text, HTML('&nbsp;'))
  #   
  #   
  # })
  
}

set_UI_clip_information <- function(query_results) {
  
  ## Get binding sites info
  encori_overlaps <- NULL
  encori_span_tag <- NULL

  ## Connect to the CLIP data database and retrieve data from the current gene
  database_path <- file.path("./database/clip_data.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  if ( DBI::dbExistsTable(conn = con, name = (query_results$gene_name %>% unique )[1]) ) {
    
    query <- paste0("SELECT * FROM '", (query_results$gene_name %>% unique )[1], "'")
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
  
  if ( (jxn_info %>% drop_na(clinvar) %>% pull(clinvar) %>% unique) == 1 ) {
    
    #print(jxn_info$clinvar)
    #print(jxn_info$clinvar_locus)
    
    if ( any(jxn_info %>% names == "clinvar_locus") ) {
      
      clinvar_hits <- (jxn_info %>% drop_na("clinvar_locus") %>% pull(clinvar_locus) %>% unique) %>%
        str_split(pattern = ",") %>%
        unlist() %>%
        unique() %>% 
        sort()
      
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
                                           paste(clinvar_hits,collapse = "; ")), HTML('&nbsp;'))))
    }
  }
  
  if ( is.null(clinvar_span_tag) ) {
    clinvar_span_tag <- div(p(strong("ClinVar variants:"), HTML('&nbsp;'), "No ClinVar data available for this junction."))
  }
  
  return(clinvar_span_tag)
}
