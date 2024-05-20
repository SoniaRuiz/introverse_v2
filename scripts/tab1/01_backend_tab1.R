###################################################
## JUNCTION QUERY - DATABASE QUERY AND UI FUNCTIONS
###################################################

query_database_nearby <- function(chr, start, end, strand,
                                  nearby_window, 
                                  databases_name, 
                                  project_id,
                                  table_name) {
  
  
  
  #chr = 10
  #start = 87880439
  #end = 87925512
  #strand = "+"
  #databases_name = "encode"
  #project_id = "TARDBP"
  #table_name = "case"
  #nearby_window = 0
  # 
  #   chr = 10
  #   start = 87932343
  #   end = 87933012
  #   strand = "+"
  #   databases_name = "encode"
  #   project_id = "TARDBP"
  #   table_name = "case"
  #   nearby_window = 100000
  
  
  
  message("'", chr, ":", start, "-", end, ":", strand, "' # ", nearby_window, " - " , databases_name, " - ", project_id, " - ", table_name)
  
  if ( any(databases_name == "all") ){
    databases_sqlite_list <- database_equivalences$sqlite_file
  } else {
    databases_sqlite_list <- database_equivalences %>% filter(key == databases_name) %>% pull(sqlite_file)
  }
  
  query_results <- map_df(databases_sqlite_list, function(database_sqlite) {
    
    # database_sqlite <- databases_sqlite_list[1]
    # database_sqlite <- databases_sqlite_list[2]
    # print(database_sqlite)
    # print(table_name)
    
    
    
    ## Connect the database
    # database_path <- file.path(database_folder, database_sqlite)
    # con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
    con <- conn_list[database_sqlite][[1]]
    
    ## Query the global database information
    
    ## Retrieve the junction type
    query <- paste0("SELECT ref_junID, ref_coordinates AS coordinates, seqnames, start, end, strand, transcript_id
                    FROM 'intron' 
                    WHERE seqnames='chr", chr, "' AND start >= ", (start %>% as.integer() - nearby_window%>% as.integer())," AND end <= ", (end %>% as.integer() + nearby_window%>% as.integer()), " AND strand='", strand, "'")
    master_intron_information <- DBI::dbGetQuery(con, query) %>% 
      arrange(start) %>%
      mutate(junction_type = "annotated_intron")
    
    
    query <- paste0("SELECT ref_junID, novel_junID, novel_coordinates AS coordinates, seqnames, start, end, strand, novel_type
                    FROM 'novel' 
                    WHERE seqnames='chr", chr, "' AND start >= ", start %>% as.integer() -nearby_window %>% as.integer()," AND end <= ", end %>% as.integer() +nearby_window %>% as.integer(), " AND strand='", strand, "'")
    master_novel_information <- DBI::dbGetQuery(con, query) %>%
      dplyr::rename(junction_type = novel_type) 
    master_novel_information <- master_novel_information %>%
      left_join(y = DBI::dbGetQuery(con, paste0("SELECT ref_junID, transcript_id FROM 'intron' WHERE ref_junID IN ('", paste(master_novel_information$ref_junID %>% unique, collapse="','"), "')")),
                by = "ref_junID")
    
    master_information <- plyr::rbind.fill(master_intron_information, master_novel_information) %>%
      arrange(start) %>%
      distinct(coordinates, .keep_all = T)
    
    
    
    ## Get the metadata
    if ( master_information %>% nrow > 0 ) {
      
      if ( project_id == "all" ) {
        
        db_misspliced_clusters <- master_information
        
      } else {
        
        if ( table_name == "all" ) {
          
          if (project_id == "all") {
            query <- paste0("SELECT * FROM 'metadata'")
          } else {
            query <- paste0("SELECT * FROM 'metadata' WHERE SRA_project='", project_id, "'")
          }
          
          db_metadata <- DBI::dbGetQuery(con, query)
          
          db_misspliced_clusters <- map_df(db_metadata$SRA_project %>% unique, function(project_id) {
            
            map_df(db_metadata %>% filter(SRA_project == project_id) %>% distinct(cluster) %>% pull, function(cluster) {
              
              message(project_id, " - ", cluster)
              
              query <- paste0("SELECT ref_junID, novel_junID FROM '", cluster, "_", project_id, "_misspliced'") 
              misspliced <- DBI::dbGetQuery(con, query) %>% 
                inner_join(y = master_information, by = c("ref_junID", "novel_junID"))
              
              query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster, "_", project_id, "_nevermisspliced'")
              nevermisspliced <- DBI::dbGetQuery(con, query) %>% 
                inner_join(y = master_information, by = "ref_junID")
              
              
              plyr::rbind.fill(misspliced, nevermisspliced) %>% distinct(ref_junID, novel_junID, .keep_all = T) %>% return()
              
            })
          })
          
        } else {
          
          query <- paste0("SELECT ref_junID, novel_junID FROM '", table_name, "_", project_id, "_misspliced'") 
          misspliced <- DBI::dbGetQuery(con, query) %>% 
            filter(ref_junID %in% master_information$ref_junID | 
                     novel_junID %in% master_information$novel_junID) 
          
          query <- paste0("SELECT DISTINCT ref_junID FROM '", table_name, "_", project_id, "_nevermisspliced'")
          nevermisspliced <- DBI::dbGetQuery(con, query) %>%  
            as_tibble() %>%
            filter(ref_junID %in% c(master_information %>% filter(junction_type == "annotated_intron") %>% distinct(ref_junID) %>% pull(ref_junID)))
          
          db_misspliced_clusters <- plyr::rbind.fill(misspliced, nevermisspliced) %>% distinct(ref_junID, novel_junID, .keep_all = T)
          
        }
        
        
        
      } 
      
      
      
      
      if ( db_misspliced_clusters %>% nrow() > 0) {
        
        db_misspliced_clusters <- rbind( db_misspliced_clusters %>%
                                           dplyr::distinct(ref_junID) %>%
                                           inner_join(y = master_information %>% filter(junction_type == "annotated_intron") %>% distinct(ref_junID, .keep_all = T),
                                                      by = "ref_junID"),
                                         db_misspliced_clusters %>%
                                           dplyr::distinct(novel_junID) %>%
                                           drop_na() %>%
                                           inner_join(y = master_information %>% filter(junction_type != "annotated_intron") %>% distinct(novel_junID, .keep_all = T),
                                                      by = "novel_junID"))
        
        ## GET GENE INFORMATION
        query <- paste0("SELECT transcript.id, transcript.transcript_id AS transcript_ENS, transcript.MANE, gene.gene_name, gene.gene_id
                        FROM 'transcript' INNER JOIN 'gene' ON transcript.gene_id=gene.id
                        WHERE transcript.id IN (", paste(db_misspliced_clusters %>% drop_na(transcript_id) %>% pull(transcript_id) %>% unique, 
                                                         collapse = ","), ")")
        db_gene <- DBI::dbGetQuery(con, query) 
        
        
        
        db_misspliced_clusters %>%
          distinct(coordinates, .keep_all = T) %>%
          inner_join(y = db_gene %>% dplyr::select(id, gene_name, gene_id),
                     by = c("transcript_id" = "id")) %>% 
          dplyr::select(-c("ref_junID","transcript_id", "novel_junID")) %>%
          return()
      }
      
    } else {
      return(NULL)
    }
    
  })
  
  # query_results %>% as.data.frame()
  
  if ( query_results %>% nrow > 0 ) {
    
    query_results %>%
      distinct(coordinates, .keep_all = T) %>%
      mutate(junction_type = stringr::str_replace_all(string = junction_type, pattern = "_", replacement = " ") %>% stringr::str_to_title()) %>%
      dplyr::rename(Coordinates = coordinates,
                    "Junction Type" = junction_type,
                    "Gene Symbol" = gene_name,
                    "Gene ID" = gene_id,
                    "Chromosome" = seqnames,
                    "Start" = start,
                    "End" = end,
                    "Strand" = strand) %>%
      return()
    
  } else {
    
    return(data.frame("Value" ="No splicing events found!"))
  }
  
}

query_database <- function(chr, start, end, strand, 
                           databases_name, 
                           project_id = NULL,
                           table_name = NULL) {
  
  # chr = 10
  # start = 87880439
  # end = 87925512
  # strand = "+"
  # databases_name = "gtex"
  # project_id = "ADIPOSE_TISSUE"
  # table_name = "Adipose - Subcutaneous"
  
  # chr = 19
  # start = 4491836
  # end = 4492014
  # strand ="+"
  # databases_name = "tcga"
  # project_id = "LAML"
  # table_name = "Primary blood derived cancer - peripheral blood"
  
  # chr = 19
  # start = 4501318
  # end = 4501913
  # strand ="+"
  # databases_name = "all"
  # project_id = "all"
  # table_name = "all"
  # chr = 19
  # start = 4491836
  # end = 4492014
  # strand ="+"
  # databases_name = "gtexv8"
  # project_id = "ADIPOSE_TISSUE"
  # table_name = "Adipose - Subcutaneous"
  #   
  # chr = 1
  # start = 154626001
  # end = 154627953
  # strand ="-"
  # databases_name = "tcga"
  # project_id = "ACC"
  # table_name = "Primary tumor"
  
  # chr = 1
  # start = 108896431
  # end = 108896864
  # strand = "+"
  
  #"chr10:87862643-87925512:+"
  # chr = 10
  # start = 87862643
  # end = 87925512
  # strand = "+"
  # databases_name = "encode"
  # project_id = "TARDBP"
  # table_name = "case"
  # 110271
  # novel_junID = 143998
  
  # chr = 1
  # start=155235312
  # end= 155235462
  # strand="-"
  
  # chr = 1
  # start = 155235308
  # end = 155235680
  # strand="-"
  # table_name = NULL
  
  # chr = 1
  # start = 154597977
  # end = 154598401
  # strand = "-"
  
  
  # chr = 1
  # start = 154597977
  # end = 154598075
  # strand = "-"
  
  
  # chr = 19
  # start=17621871
  # end=17623117
  # strand="-"
  
  # chr = 16
  # start = 86478133
  # end = 86508654
  # strand = "-"
  
  # databases_name = "all"
  
  # databases_name = "tcga"
  # project_id = "COAD"
  # table_name = "all"
  # table_name = "Additional metastatic"
  
  # databases_name = "GTEx v8"
  # databases_name = "ENCODE shRNA"
  # databases_name = "Alzheimer's Disease/Control"
  # databases_name = "Brainseq Squizo/Control"
  # databases_name = "Parkinson's Disease/Control"
  # databases_name = "LRRK2 G2019S mutation"
  
  # print(databases_name)
  
  if ( any(databases_name == "all") ){
    databases_sqlite_list <- database_equivalences$sqlite_file
    project_id <- NULL
    table_name <- NULL
  } else {
    databases_sqlite_list <- database_equivalences %>% filter(key == databases_name) %>% pull(sqlite_file)
  }
  
  
  
  query_results <- map_df(databases_sqlite_list, function(database_sqlite) {
    
    # database_sqlite <- databases_sqlite_list[1]
    # database_sqlite <- databases_sqlite_list[2]
    # print(database_sqlite)
    # print(table_name)
    
    ## Connect the database
    #database_path <- file.path(database_folder, database_sqlite)
    #con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
    con <- conn_list[database_sqlite][[1]]
    
    ## Get the metadata
    if ( is.null(project_id) || project_id == "all" ) {
      query <- paste0("SELECT * FROM 'metadata'")
    } else {
      if ( is.null(table_name) || table_name == "all") {
        query <- paste0("SELECT * FROM 'metadata' WHERE SRA_project='", project_id, "'")
      } else {
        query <- paste0("SELECT * FROM 'metadata' WHERE SRA_project='", project_id, "' AND cluster='",table_name,"'")
      }
    }
    db_metadata <- DBI::dbGetQuery(con, query)
    
    
    ## Set the junID coordinates
    jxn_of_interest_coordinates <- paste0("chr", chr, ":", start, "-", end, ":", strand)
    print(jxn_of_interest_coordinates)
    
    ## Query the global database information
    
    
    ## Retrieve the junction type
    query <- paste0("SELECT * 
                    FROM 'intron' 
                    WHERE ref_coordinates='", jxn_of_interest_coordinates, "'")
    master_information <- DBI::dbGetQuery(con, query) 
    
    ## It is not an annotated intron. We check within the master novel table
    if ( (master_information %>% nrow) == 0 ) {
      
      query <- paste0("SELECT * 
                      FROM 'novel' 
                      WHERE novel_coordinates='", jxn_of_interest_coordinates, "'")
      master_information <- DBI::dbGetQuery(con, query)
      
      ## It is a novel junction
      if ( (master_information %>% nrow) > 0 ) {
        
        junction_type = master_information$novel_type %>% unique()
        
        ## We add information about the reference intron
        query <- paste0("SELECT ref_junID, ref_coordinates, transcript_id_list 
                        FROM 'intron' 
                        WHERE ref_junID='", master_information$ref_junID, "'")
        intron_coordinates <- DBI::dbGetQuery(con, query) 
        
        master_information <- master_information  %>% 
          dplyr::rename(length = novel_length,
                        mes5ss = novel_mes5ss,
                        mes3ss = novel_mes3ss,
                        donor_sequence = novel_donor_sequence,
                        acceptor_sequence = novel_acceptor_sequence)%>%
          dplyr::mutate(junID = "novel_junID") %>%
          left_join(intron_coordinates,
                    by = "ref_junID")
      } 
      
    } else {
      
      junction_type = "annotated intron"
      
      ## We add info about the mis-spliced junctions attached to the annotated intron of interest
      query <- paste0("SELECT ref_junID, novel_junID, novel_coordinates, novel_type 
                      FROM 'novel' 
                      WHERE ref_junID='", master_information$ref_junID, "'")
      novel_coordinates <- DBI::dbGetQuery(con, query) 
      
      master_information <- master_information %>% 
        dplyr::rename(length = ref_length,
                      mes5ss = ref_mes5ss,
                      mes3ss = ref_mes3ss,
                      donor_sequence = ref_donor_sequence,
                      acceptor_sequence = ref_acceptor_sequence) %>%
        dplyr::select(-transcript_id) %>%
        left_join(novel_coordinates,
                  by = "ref_junID") %>%
        dplyr::mutate(junID = "ref_junID") 
    }
    
    if ( (master_information %>% distinct(ref_junID) %>% nrow) == 1) {
      
      ## Loop per each database table to get the actual splicing info 
      ## for the junction of interest in the current cluster
      local_query_result <- map_df(db_metadata$cluster %>% unique(), function(cluster_id) {
        
        # cluster_id <- (db_metadata$cluster %>% unique())[1]
        
        message(cluster_id, " - retrieving splicing information...")
        
        project_ids <- db_metadata %>%
          filter(cluster == cluster_id) %>%
          distinct(SRA_project) %>%
          pull()
        
        map_df(project_ids, function(project_id) {
          
          message(project_id, "...")
          
          # project_id <- project_ids[1]
          # project_id <- project_ids[19]
          
          ## We get the actual splicing info for the junction of interest in the current cluster
          if ( junction_type == "annotated intron" ) {
            
            query <-  paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced' 
                             WHERE ref_junID='", master_information$ref_junID %>% unique, "'")
            introns <- dbGetQuery(con, query) 
            query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_nevermisspliced' 
                            WHERE ref_junID='", master_information$ref_junID %>% unique, "'" )
            introns <- plyr::rbind.fill(introns, dbGetQuery(con, query)) 
            
          } else {
            
            query <-  paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced' 
                             WHERE novel_junID='", master_information$novel_junID, "'")
            introns <- dbGetQuery(con, query) 
            
          }
          
          if ( introns %>% nrow > 0 ) {
            ## We join the mis-splicing info from the current cluster and the info from the master table
            introns %>% 
              mutate(cluster = cluster_id, 
                     project = project_id,
                     coordinates = jxn_of_interest_coordinates,
                     junID_type = junction_type,
                     database = database_equivalences %>% filter(sqlite_file == database_sqlite) %>% pull(name),
                     cluster_index = database_equivalences %>% filter(sqlite_file == database_sqlite) %>% pull(index),
                     cluster_type = database_equivalences %>% filter(sqlite_file == database_sqlite) %>% pull(project_types) ) %>%
              return()
          } else {
            return(NULL)
          }
          
        })
        
      })
      
      if ( local_query_result %>% nrow() > 0) {
        
        ## We add the master information about the genes and transcripts
        query <- paste0("SELECT transcript.id, transcript.transcript_id AS transcript_ENS, transcript.MANE, transcript.transcript_biotype, gene.gene_name 
                        FROM 'transcript' INNER JOIN 'gene' ON transcript.gene_id=gene.id
                        WHERE transcript.id IN (", master_information$transcript_id_list %>% unique, ")")
        db_gene <- DBI::dbGetQuery(con, query)
        
        if ( junction_type == "annotated intron" && all(local_query_result$ref_type %>% unique == "never") ) {
          local_query_result <- local_query_result %>% dplyr::select(-novel_junID) %>%
            inner_join(y = master_information,
                       by = c("ref_junID" = "ref_junID"))
        } else {
          local_query_result <- local_query_result %>%
            left_join(y = master_information,
                      by = c("ref_junID" = "ref_junID",
                             "novel_junID" = "novel_junID")) 
        }
        local_query_result %>%
          full_join(y = db_gene,
                    by = c("transcript_id" = "id") ) %>% 
          return()
      } else {
        print("No results found for the criteria selected!")
        return(NULL)
      } 
      
      } else {
        print("No results found for the criteria selected!")
        return(NULL)
        
      }
    
    # } else {
    #   print("No results found!")
    #   return(NULL)
    # }  
    
    
  })
  
  #DBI::dbDisconnect(con)
  query_results %>%
    as_tibble() %>%
    return()
  
  
  
  }


setup_UI_details_section <- function(query_results, back_button = NULL) {
  
  # message("Back button...", back_button %>% as.integer())
  
  v <- list()
  i <- 1
  
  common_jxn_info <- query_results %>% distinct(coordinates, transcript_id, .keep_all = T) %>% drop_na(ref_junID)
  
  
  # Set UI data
  
  distance_and_global_info <- set_UI_distance_information(jxn_info = common_jxn_info)
  
  transcript_span_tag <- set_UI_transcript_information(query_results)
  
  encori_span_tag <- set_UI_clip_information(query_results)
  
  clinvar_span_tag <- set_UI_clinvar_information(jxn_info = common_jxn_info)
  
  
  if ( !is.null(back_button) && back_button != "" && back_button %>% as.integer() > 0) {
    back_button = div(tags$a(href = "#", role="button",  onclick = "var value=$('#nearby_helper_tab1').val();Shiny.setInputValue('database_nearby_tab1', value );Shiny.setInputValue('nearby_helper_tab1', '');Shiny.onInputChange('acceptButton_tab1', 1 ,{priority:'event'});", "< Back to results"))
  } else {
    back_button = span("")
  }
  
  ## Return UI data
  
  v[[i]] <- div(
    back_button,
    ## Set hidden info
    p(id = "phastCons17_5ss", style = "display:none;", common_jxn_info %>% drop_na(mean_phastCons17way5ss_100) %>% pull(mean_phastCons17way5ss_100) %>% unique),
    p(id = "phastCons17_3ss", style = "display:none;", common_jxn_info %>% drop_na(mean_phastCons17way3ss_100) %>% pull(mean_phastCons17way3ss_100) %>% unique),
    p(id = "CDTS_5ss", style = "display:none;", common_jxn_info %>% drop_na(mean_CDTS5ss_100) %>% pull(mean_CDTS5ss_100) %>% unique),
    p(id = "CDTS_3ss", style = "display:none;", common_jxn_info %>% drop_na(mean_CDTS3ss_100) %>% pull(mean_CDTS3ss_100) %>% unique),
    p(id = "MES_5ss", style = "display:none;", common_jxn_info %>% drop_na(mes5ss) %>% pull(mes5ss) %>% unique),
    p(id = "MES_3ss", style = "display:none;", common_jxn_info %>% drop_na(mes3ss) %>% pull(mes3ss) %>% unique),
    p(id = "donor_sequence", style = "display:none;", common_jxn_info %>% drop_na(donor_sequence) %>% pull(donor_sequence) %>% unique),
    p(id = "acceptor_sequence", style = "display:none;", common_jxn_info %>% drop_na(acceptor_sequence) %>% pull(acceptor_sequence) %>% unique),
    

    
    p(id = "CLNSIG_list", style = "display:none;", common_jxn_info %>% drop_na(CLNSIG_list) %>% pull(CLNSIG_list)),
    p(id = "CLNVC_list", style = "display:none;", common_jxn_info %>% drop_na(CLNVC_list) %>% pull(CLNVC_list)),
    p(id = "MC_list", style = "display:none;", common_jxn_info %>% drop_na(MC_list) %>% pull(MC_list)),
    
    
    ## Set visible info
    h3("Junction Details"),
    hr(),
    p(strong("Coordinates:"), span(common_jxn_info$coordinates %>% unique, id = "span_junID")),
    p(strong("Length (in bp):"), span(common_jxn_info %>% drop_na(length) %>% pull(length) %>% unique)),
    p(strong("Junction type:", 
             span(str_to_title(string = str_replace(string = common_jxn_info$junID_type %>% unique, pattern = "_", replacement = " ")), 
                  style = paste0("color:", distance_and_global_info$color[1] %>% unique), 
                  id="span_junType")),
      " (Ensembl v111)."),
    p(strong(distance_and_global_info$distance_label[2]), distance_and_global_info$distance_info[2]), ## Ref junction
    p(strong(distance_and_global_info$distance_label[1], 
             span(distance_and_global_info$distance_info[1], 
                  style = paste0("color:", distance_and_global_info$color[2] %>% unique)))), ## Distance
    p(strong("Gene:"), span(paste(common_jxn_info$gene_name %>% unique, collapse = " "), id = "span_geneName")),
    p(strong("Transcript visualisation:"), HTML('&nbsp;'), transcript_span_tag, id="transcript_list"),
    encori_span_tag,
    clinvar_span_tag,
    p(tags$a(href = "#", onclick = "
             
             Shiny.setInputValue('phastCons17_5ss_tab1', $('#phastCons17_5ss')[0].innerText);
             Shiny.setInputValue('phastCons17_3ss_tab1', $('#phastCons17_3ss')[0].innerText);
             Shiny.setInputValue('CDTS_5ss_tab1', $('#CDTS_5ss')[0].innerText);
             Shiny.setInputValue('CDTS_3ss_tab1', $('#CDTS_3ss')[0].innerText);
             Shiny.setInputValue('MES_5ss_tab1', $('#MES_5ss')[0].innerText);
             Shiny.setInputValue('MES_3ss_tab1', $('#MES_3ss')[0].innerText);
             Shiny.setInputValue('donor_sequence_tab1', $('#donor_sequence')[0].innerText);
             Shiny.setInputValue('acceptor_sequence_tab1', $('#acceptor_sequence')[0].innerText);


             
             $('#modalVisualiseMORE_tab1').modal('show');", 
             "More info...")),
    class = "well"
  )
  
  return(v)
}


setup_UI_results_section <- function(query_results) {
  
  
  databases_queried <- query_results %>%
    distinct(database) %>%
    drop_na() %>%
    pull() 
  
  result_tabs <- list()
  
  for ( j in 1:(length(databases_queried)) ) {
    
    # j = 1
    # j = 2
    
    local_query_results <- query_results %>% 
      filter(database == databases_queried[j]) %>% 
      group_by(project, cluster) %>%
      distinct( ref_junID, novel_junID, .keep_all = T) %>%
      ungroup()
    
    
    current_db_general_info <- database_equivalences %>% filter(name == local_query_results$database %>% unique)
    
    # print(current_db_general_info)
    
    #######################################################################################
    # 1. Get avg number of individuals and reads corresponding to the junx of interest
    #######################################################################################
    
    local_query_results_tidy_stats <- local_query_results %>% 
      group_by(project, cluster) %>%
      distinct(coordinates, .keep_all = T) %>%
      ungroup()
    
    if ( local_query_results_tidy_stats$junID_type %>% unique == "annotated intron" ) {
      n_individuals = local_query_results_tidy_stats$ref_n_individuals %>% sum 
      n_reads = local_query_results_tidy_stats$ref_sum_counts %>% mean %>% round(digits = 2)
    } else {
      n_individuals = local_query_results_tidy_stats$novel_n_individuals %>% sum 
      n_reads = local_query_results_tidy_stats$novel_sum_counts %>% mean %>% round(digits = 2)
    }
    
    
    #######################################################################################
    # 2. Get message of the number clusters and projects were the junction has been per database
    #######################################################################################
    
    database_metadata <- get_database_metadata(database_sqlite_list = current_db_general_info$sqlite_file)
    
    per_individuals_message = paste0(n_individuals, " samples (", 
                                     round((n_individuals * 100)/
                                             (database_metadata[[ 1 ]] %>% 
                                                group_by(SRA_project) %>% 
                                                distinct(sample_id) %>% 
                                                nrow()), digits=2),
                                     "% of ", 
                                     (database_metadata[[ 1 ]] %>% 
                                        group_by(SRA_project) %>% 
                                        distinct(sample_id) %>% 
                                        nrow())," total samples studied in this database)")
    
    
    #######################################################################################
    ## Tidy database results prior returning to UI
    #######################################################################################
    
    ## 1. Remove columns depending on the type of the junction queried
    if ( local_query_results$junID_type %>% unique == "annotated intron" ) {
      
      local_query_results <- local_query_results %>%
        dplyr::relocate(dplyr::any_of(x = c("project",
                                            "cluster",
                                            "ref_coordinates", 
                                            "ref_n_individuals",
                                            "ref_sum_counts",
                                            "MSR_D","MSR_A",
                                            "novel_coordinates", "novel_type")))
      
      ## If it is a NEVER-MISSPLICED annotated intron, we remove data corresponding to the novel junction
      if ( all(local_query_results$ref_type %>% unique == "never") ) {
        
        local_query_results <- local_query_results %>%
          dplyr::select(-c(novel_coordinates, 
                           novel_n_individuals,
                           novel_sum_counts,
                           novel_type)) %>%
          distinct(database,project,cluster,ref_junID,.keep_all = T)
        
      } else {
        local_query_results <- local_query_results %>%
          mutate(ref_type = paste0(ref_type, " splice sites")) 
      }
      
    } else {
      ## If it is a novel junction, we remove data corresponding to the annotated reference intron
      local_query_results <- local_query_results %>%
        dplyr::relocate("project","cluster","novel_coordinates","novel_n_individuals", "novel_sum_counts", "novel_type", "ref_coordinates",
                        "ref_n_individuals","ref_sum_counts","ref_type") %>%
        dplyr::select(-any_of(x = c("ref_donor_sequence",
                                    "ref_acceptor_sequence","ref_mes5ss","ref_mes3ss","MSR_D","MSR_A")))
    }
    
    
    #######################################################################################
    ## 2. Rename junctions
    #######################################################################################
    
    rename_column_lookup <- c("Intron Coordinates" = "ref_coordinates",
                              "Novel Junction Coordinates" = "novel_coordinates",
                              "Intron seen in no. samples" = "ref_n_individuals",
                              "Intron Total Reads" = "ref_sum_counts",
                              "Mis-spliced At" = "ref_type",
                              "MSR Donor" = "MSR_D",
                              "MSR Acceptor" = "MSR_A",
                              "Novel Junction Num. Samples" = "novel_n_individuals",
                              "Novel Junction Total Reads" = "novel_sum_counts",
                              "Junction Type" = "novel_type",
                              "Sample Group" = "cluster",
                              "Gene Name" = "gene_name",
                              "Gene TPM" = "gene_tpm",
                              "Project" = "project",
                              "Database" = "database")
    
    final_query_results <- local_query_results %>%
      mutate(cluster = ifelse(cluster == "case", "shRNA knockdown", cluster)) %>%
      dplyr::relocate(gene_tpm, .before = gene_name) %>%
      
      mutate_at(vars(one_of('novel_type')), ~str_to_sentence(str_replace(string = .x, pattern = "_", replacement = " "))) %>%
      
      mutate(project = str_replace(string = project, pattern = "_", replacement = " "),
             across(where(is.numeric), round, digits = 5)) %>%
      
      dplyr::select(-c(dplyr::any_of(x = c("ref_junID","novel_junID.x", "novel_junID.y", "transcript_id", "transcript_id_list", "novel_junID",
                                           "seqnames", "start", "end", "strand", "distance", "junID", "coordinates", "cluster_index", "misspliced",
                                           "u2_intron", "protein_coding", "cluster_type", "junID_type","transcript_ENS",
                                           "novel_donor_sequence","novel_acceptor_sequence",
                                           "mean_phastCons17way5ss_100","mean_phastCons17way3ss_100",
                                           "mean_CDTS5ss_100","mean_CDTS3ss_100", "clinvar","ref_type", "clinvar_locus",
                                           "donor_sequence", "acceptor_sequence", "length","database","gene_name","MANE",
                                           "transcript_biotype", 
                                           "CLNSIG_list",
                                           "CLNVC_list",
                                           "MC_list",
                                           "mes5ss","mes3ss"))))  %>% 
      
      dplyr::rename(dplyr::any_of(rename_column_lookup)) 
    
    
    
    #######################################################################################
    ## 4. Return data
    #######################################################################################
    
    result_tabs[[j]] <-  tabPanel(title = local_query_results$database %>% unique(),
                                  id = j,
                                  div(
                                    ## Set visible info
                                    h4(id="h4_databaseID", strong(local_query_results$database %>% unique())),
                                    p(strong("Found in: "), span(per_individuals_message,  style = "color:blue")),
                                    p(strong("Avg. reads: "), span(n_reads, style = "color:blue")),
                                    DT::datatable( data =  final_query_results,
                                                   rownames = F,
                                                   extensions = c('Buttons', 'Responsive'),
                                                   options = list(
                                                     lengthMenu = c(5, 10, 20, 40, 80, 100),
                                                     pageLength = 15,
                                                     scrollCollapse=TRUE,
                                                     dom = '<"top"fB>t<"bottom"lpi><"clear">',
                                                     buttons = c('copy', 'csv', 'excel')
                                                   ))
                                  ),
                                  hr()   
    )
  }
  
  return(result_tabs)
}
