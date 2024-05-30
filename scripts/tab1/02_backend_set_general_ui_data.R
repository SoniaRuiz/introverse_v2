


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
      if ( any(jxn_info$ref_type %>% unique == "both" ) || 
           (any(jxn_info$ref_type %>% unique == "donor" ) && 
            any(jxn_info$ref_type %>% unique == "acceptor" )) ) {
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
  
  
  if ( str_detect(query_results %>% drop_na(junID_type) %>% pull(junID_type) %>% unique, pattern = "novel") ) {
    junID = query_results %>% drop_na(novel_coordinates) %>% pull(novel_coordinates) %>% unique
  } else {
    junID = query_results %>% drop_na(ref_coordinates) %>% pull(ref_coordinates) %>% unique
  }

  
  query_results_transcript <- get_transcript_to_plot(junID = junID,
                                                     geneName = query_results %>% drop_na(gene_name) %>% pull(gene_name) %>% unique,
                                                     jxn.type = query_results %>% drop_na(junID_type) %>% pull(junID_type) %>% unique,
                                                     multiple = T)
  

  transcript_span_tag = NULL

  
  ## Detect if the junction can be plotted in the MANE transcript and in any of the other protein-coding transcripts
  ## It might happen that the junction is part of a transcript with a non-protein-coding biotype, hence cannot be plotted.
  
  if ( query_results_transcript$transcript %>% drop_na(tag) %>% nrow() == 1 ) {

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
      $('#modalVisualiseTranscript_tab1').modal('show');", 
               query_results_transcript$transcript %>% drop_na(tag) %>% pull(transcript_id) %>% unique), " (MANE)")

  } 
  
    
  if ( query_results_transcript$utr %>% filter(tag %>% is.na()) %>% distinct(transcript_id) %>% nrow > 0 ) {
    
    if ( !is.null(transcript_span_tag) ) { 
      transcript_span_tag <- span(transcript_span_tag, "and")
    }
    
    transcript_span_tag <- span(transcript_span_tag, 
                                span(tags$a(href = "#", 
                                onclick = "var ID_jxn= $('#span_junID')[0].innerText; 
                                Shiny.setInputValue('junID_tab1', ID_jxn);
                                
                                var ID_database= $('#h4_databaseID')[0].innerText; 
                                Shiny.setInputValue('databaseName_tab1', ID_database);
                                
                                //var transcript_ID = $('#atag_transcriptID')[0].innerText;
                                //Shiny.setInputValue('transcriptENS_tab1', transcript_ID );
  
                               var junType= $('#span_junType')[0].innerText; 
                               Shiny.setInputValue('junType_tab1', junType);
  
                               var gene_Name= $('#span_geneName')[0].innerText; 
                               Shiny.setInputValue('geneName_tab1', gene_Name);
  
                               Shiny.setInputValue('transcriptMANE_tab1', 0);
                               $('#modalVisualiseTranscript_tab1').modal('show');", 
                                query_results_transcript$utr %>% filter(tag %>% is.na()) %>% distinct(transcript_id) %>% nrow(), 
                                paste(str_replace_all(string = query_results_transcript$utr %>% filter(tag %>% is.na()) %>% distinct(transcript_id, .keep_all = T) %>% pull(transcript_biotype) %>% unique,
                                                      pattern = "_", replacement = " "), collapse = "' and '"), 
                                "transcripts"
                               )))
    
  }
  
  if ( is.null(transcript_span_tag) ) {
    transcript_span_tag <- div(p(strong("Transcript visualisation:"), HTML('&nbsp;'), "There are no protein-coding transcripts available for this gene (Ensembl 111)."))
  } else {
    transcript_span_tag <- div(p(strong("Transcript visualisation:"), HTML('&nbsp;'), transcript_span_tag))
  }
  
  return(transcript_span_tag)
  

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
                                                       dplyr::select(seqnames, start, end, strand ) %>% 
                                                       drop_na() %>% 
                                                       distinct(.keep_all = T)  %>%
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
                                       //var transcript_Name= $('#transcript_list span a')[0].innerText; 
                                       //Shiny.setInputValue('transcriptENS_tab1', transcript_Name );
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
    ## Get clinvar overlaps
    
    ## Connect to the CLIP data database and retrieve data from the current gene
    clinvar_data <- get_database_clinvar(jxn_info) 
    
    if ( clinvar_data %>% nrow() > 0 ) {
      
      
      clinvar_span_tag <-div(p(strong("ClinVar variants:"), HTML('&nbsp;'),
                               span(tags$a(href = "#", 
                                           onclick = "var ID_jxn= $('#span_junID')[0].innerText; 
                                           Shiny.setInputValue('junID_tab1', ID_jxn);
                                           
                                           var ID_database= $('#h4_databaseID')[0].innerText; 
                                           Shiny.setInputValue('databaseName_tab1', ID_database);

                                           //Shiny.setInputValue('clinvarlocus_tab1', $('#clinvar_list')[0].innerText );

                                           var junType= $('#span_junType')[0].innerText; 
                                           Shiny.setInputValue('junType_tab1', junType);

                                           var gene_Name= $('#span_geneName')[0].innerText; 
                                           Shiny.setInputValue('geneName_tab1', gene_Name);
                                           
                                           $('#modalVisualiseClinVar_tab1').modal('show');",
                                           
                                           clinvar_data$ID  %>% unique() %>% length(), " splicing variants"
                                           
                                           ), HTML('&nbsp;'))))
    }
  }
  
  if ( is.null(clinvar_span_tag) ) {
    clinvar_span_tag <- div(p(strong("ClinVar variants:"), HTML('&nbsp;'), "No ClinVar data available for this junction."))
  }
  
  return(clinvar_span_tag)
}
