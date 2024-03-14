library(tidyverse)
library(DBI)
library(ggplot2)
library(ggtranscript)

#setwd()

database_folder <- base::file.path("database")

database_equivalences <- data.frame(index = c(1,2,3,4,5),
                                    key = c("gtex_v8", "gtex_v8_age", "encode_shrna", "ad_control", "brainseq"),
                                    name = c("GTEx v8", "GTEx v8 Age", "ENCODE shRNA", "AD/Control", "Brainseq Squizophrenia/Control" ),
                                    sqlite_file = c("splicing_1read.sqlite","splicing_1read_age.sqlite","ENCODE_SR_1read.sqlite","SRP100948_1read.sqlite","brainseq_1read.sqlite"),
                                    project_types = c("body sites", "body sites","RBPs","Alzheimer's Disease/Control","Brainseq Squizophrenia/Control"))
                                                         


#######################################
## FUNCTIONS
#######################################

get_database_metadata <- function(database_list) {
  
  purrr::map(database_list, function(database_name) {
    
    # database_name <- database_list[2]
    database_path <- file.path(database_folder, database_name)
    
    con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
    query <- paste0("SELECT * FROM 'metadata'")
    metadata <- DBI::dbGetQuery(con, query)
    
    DBI::dbDisconnect(conn = con)
    
    metadata %>%
      return()
    
  })
  
}


query_database <- function(chr, start, end, strand, 
                           databases_name, 
                           table_name = NULL) {
  

  # chr = 19
  # start = 44905842
  # end = 44906601
  # strand = "+"
  
  # chr = 10
  # start = 87862643
  # end = 87864453
  # strand = "+"
  
  # chr = 1
  # start=155235312
  # end= 155235462
  # strand="-"

  chr = 1
  start = 155235308
  end = 155235680
  strand="-"
  table_name = NULL
   
  # chr = 1
  # start = 154597977
  # end = 154598401
  # strand = "-"

  # chr = 19
  # start=17621871
  # end=17623541
  # strand="-"
  
  # chr = 16
  # start = 86478133
  # end = 86508654
  # strand = "-"
  
  # databases_name = "All"
  # databases_name = "GTEx v8"
  # databases_name = "ENCODE shRNA"
  # databases_name = "AD/Control"
  # databases_name = "Brainseq Squizo/Control"
  
  # print(databases_name)
  
  if ( any(databases_name == "All") ){
    databases_sqlite_list <- database_equivalences$sqlite_file
  } else {
    databases_sqlite_list <- database_equivalences %>% filter(name == databases_name) %>% pull(sqlite_file)
  }
  
  query_results <- map_df(databases_sqlite_list, function(database_sqlite) {
   
    # database_sqlite <- databases_sqlite_list[1]
    # database_sqlite <- databases_sqlite_list[2]
    
    ## Connect the database
    database_path <- file.path(database_folder, database_sqlite)
    
    con <- DBI::dbConnect(RSQLite::SQLite(), database_path)

    ## Get the metadata
    query <- paste0("SELECT * FROM 'metadata'")
    db_metadata <- DBI::dbGetQuery(con, query)
    
    ## Set the junID coordinates
    jxn_of_interest_coordinates <- paste0("chr", chr, ":", start, "-", end, ":", strand)
    
    ## If the table has not been specified, we query the database globally
    if ( is.null(table_name) ) {
      
      ## Retrieve the junction type
      query <- paste0("SELECT * FROM 'intron' WHERE ref_coordinates='", jxn_of_interest_coordinates, "'")
      master_information <- DBI::dbGetQuery(con, query) %>% as_tibble()
      
      ## It is not an annotated intron. We check within the master novel table
      if ( (master_information %>% nrow) == 0 ) {
        
        query <- paste0("SELECT * FROM 'novel' WHERE novel_coordinates='", jxn_of_interest_coordinates, "'")
        master_information <- DBI::dbGetQuery(con, query)
        
        ## It is a novel junction
        if ( (master_information %>% nrow) > 0 ) {
          
          junction_type = master_information$novel_type %>% unique()
          
          ## We add information about the reference intron
          query <- paste0("SELECT ref_junID, ref_coordinates, transcript_id_list FROM 'intron' 
                          WHERE ref_junID='", master_information$ref_junID, "'")
          intron_coordinates <- DBI::dbGetQuery(con, query) %>% as_tibble()
          
          master_information <- master_information %>%
            dplyr::mutate(junID = "novel_junID") %>%
            left_join(intron_coordinates,
                      by = "ref_junID")
        } 
        
      } else {
        
        junction_type = "annotated intron"
        
        ## We add info about the mis-spliced junctions attached to the annotated intron of interest
        query <- paste0("SELECT ref_junID, novel_junID, novel_coordinates, novel_type FROM 'novel' 
                        WHERE ref_junID='", master_information$ref_junID, "'")
        novel_coordinates <- DBI::dbGetQuery(con, query) %>% as_tibble()
        
        master_information <- master_information %>% 
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
            
            # project_id <- project_ids[1]
            # project_id <- project_ids[19]
            
            ## We get the actual splicing info for the junction of interest in the current cluster
            if (junction_type == "annotated intron") {
              
              query <-  paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced' 
                               WHERE ref_junID='", master_information$ref_junID %>% unique, "'")
              introns <- dbGetQuery(con, query) %>% as_tibble()
              query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_nevermisspliced' 
                              WHERE ref_junID='", master_information$ref_junID %>% unique, "'" )
              introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble()) 
              
            } else {
              
              query <-  paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced' 
                               WHERE novel_junID='", master_information$novel_junID, "'")
              introns <- dbGetQuery(con, query) %>% as_tibble()
              
            }
            
            if (introns %>% nrow > 0 ) {
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
          query <- paste0("SELECT transcript.id, transcript.transcript_id AS transcript_ENS, gene.gene_name 
                          FROM 'transcript' INNER JOIN 'gene' ON transcript.gene_id=gene.id
                          WHERE transcript.id IN (", master_information$transcript_id_list %>% unique, ")")
          db_gene <- DBI::dbGetQuery(con, query)
          
          local_query_result %>%
            inner_join(y = master_information,
                       by = c("ref_junID" = "ref_junID",
                              "novel_junID" = "novel_junID")) %>%
            full_join(y = db_gene,
                       by = c("transcript_id" = "id") ) %>% 
            return()
        } else {
          return(NULL)
        } 
        
      } else {
        return(NULL)
        
      }
      
    } else {
      return(NULL)
    }  
    
    
  })
  
  #DBI::dbDisconnect(con)
  query_results %>%
    as_tibble() %>%
    return()

  
  
}

visualise_transcript <- function(junID,
                                 transcript_ENS,
                                 database_name) {
  
  # junID = "chr19:44905842-44906601:+"
  # transcript_ENS = "ENST00000446996"
  
  # junID = "chr1:154585866-154586180:-"
  # transcript_ENS = "ENST00000681786"
  # database_name = "GTEx v8"
  
  message(junID, " - ", transcript_ENS, " - ", database_name)
  
  
  ##############################
  ## GET THE TABLE NAME
  ##############################
  
  chr_junc <- junID %>%
    str_sub(start = 1,
            end = str_locate_all(string = junID, pattern = ":")[[1]][1,2]-1)
  start_junc <- junID %>%
    str_sub(start = str_locate_all(string = junID, pattern = ":")[[1]][1,2]+1,
            end = str_locate_all(string = junID, pattern = "-")[[1]][1,2]-1)
  end_junc <- junID %>%
    str_sub(start = str_locate_all(string = junID, pattern = "-")[[1]][1,2]+1,
            end = str_locate_all(string = junID, pattern = ":")[[1]][2,2]-1)
  strand_junc <- junID %>%
    str_sub(start = str_locate_all(string = junID, pattern = ":")[[1]][2,2]+1,
            end = junID %>% stringr::str_count())

  
  df_intron <- data.frame(seqnames = chr_junc %>% as.factor(),
                          start = start_junc %>% as.integer(),
                          end = end_junc %>% as.integer(), 
                          strand = strand_junc %>% as.factor(),
                          ID = junID) %>% as_tibble()

  
  # databaseID <- database_list_names[database_name] %>% unname()
  # 
  # 
  # ## Connect the database
  # database_path <- file.path(database_folder, databaseID)
  # con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  # 
  # ## Get the metadata
  # query <- paste0("SELECT * 
  #                 FROM 'metadata'
  #                 WHERE cluster ='",cluster_id,"' AND SRA_project='",project_id,"'")
  # db_metadata <- DBI::dbGetQuery(con, query)
  # 
  # 
  # ###################################
  # ## QUERY THE DATABASE FOR THE CURRENT INTRON/NOVEL JUNCTION
  # ###################################
  # 
  # 
  # 
  # sql_statement <- paste0("SELECT *
  #                         FROM 'intron' 
  #                         INNER JOIN 'transcript' ON transcript.id=intron.transcript_id
  #                         WHERE intron.ref_coordinates == '", junID, "' 
  #                         AND transcript.transcript_id='", transcript_ENS, "'")
  # 
  # df_intron <- dbGetQuery(con, sql_statement)
  #   
  # if ( df_intron %>% nrow == 0 ) ) {
  #     
  # 
  # } else {
  #   df_intron %>%
  #     dplyr::select(start, end, 
  #                   coordinates = "ref_coordinates",
  #                   n_individuals = "ref_n_individuals",
  #                   sum_counts = "ref_sum_counts") %>%
  #     distinct(coordinates,.keep_all = T)
  # }
  # else if ( str_detect(jxn_type, pattern = "annotated" ) ) {
  # 
  #   sql_statement <- paste0("SELECT *
  #                           FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue
  #                           INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
  #                           INNER JOIN 'transcript' ON transcript.id=tissue.transcript_id
  #                           WHERE tissue.ref_junID == '", junID, "' AND transcript.transcript_id='", transcript_id, "'")
  # 
  #   df_intron <- dbGetQuery(con, sql_statement)
  #   sql_statement <- paste0("SELECT *
  #                           FROM '", cluster_id, "_", project_id, "_nevermisspliced' AS tissue
  #                           INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
  #                           INNER JOIN 'transcript' ON transcript.id=tissue.transcript_id
  #                           WHERE tissue.ref_junID == '", junID, "'")
  #   df_intron <- plyr::rbind.fill(df_intron, dbGetQuery(con, sql_statement)) %>%
  #     dplyr::select(start, end, 
  #                   coordinates = "ref_coordinates",
  #                   n_individuals = "ref_n_individuals",
  #                   sum_counts = "ref_sum_counts") %>%
  #     distinct(coordinates,.keep_all = T)
  # }
    
  #######################################
  ## GET TRANSCRIPT DATA
  #######################################
  
  # hg38 <- rtracklayer::import(con = "~/PROJECTS/splicing-accuracy-manuscript/dependencies/Homo_sapiens.GRCh38.105.chr.gtf")
  # hg38_transcripts <- hg38 %>%
  #   as_tibble() %>%
  #   filter(type %in% c("transcript", "exon", "CDS", "five_prime_utr", "three_prime_utr")) %>%
  #   dplyr::select(seqnames,start,end, width, type,strand, transcript_id, transcript_biotype )
  # saveRDS(object = hg38_transcripts, file = "./dependencies/hg38_transcripts.105.rds" )
  
  hg38_transcript_to_plot <- hg38_transcripts %>%
    filter(transcript_id == transcript_ENS)
  
  
  cds_to_plot <- hg38_transcript_to_plot %>% dplyr::filter(type == "CDS")
  utr_to_plot <- hg38_transcript_to_plot %>% dplyr::filter(str_detect(type,pattern ="utr"))
  exons_to_plot <- hg38_transcript_to_plot %>% dplyr::filter(type == "exon")
  
  if (cds_to_plot %>% nrow() == 0 || utr_to_plot %>% nrow() == 0) {
    
    ggplot() +
      theme_void(base_size = 14) +
      geom_text(aes(0,0,
                    label=paste0("Sorry, it has not been possible to plot the selected junction.\nThe transcript might not have CDS or UTR regions defined.\nTranscript biotype (Ensembl v105): '", 
                                 str_replace(string = hg38_transcript_to_plot$transcript_biotype %>% unique,
                                             pattern = "_", replacement = " "), "'."))) %>%
      return()
    
  } else {
    
    #############################
    ## FIND THE REF INTRON 
    ## PRIOR THE ZOOM
    #############################
    
    all_introns <- ggtranscript::to_intron(exons_to_plot)
    
    ref_intron <- rbind(all_introns[which.min(abs(all_introns$start - df_intron$start)),],
                        all_introns[which.min(abs(all_introns$end - df_intron$end)),])
    
    ref_exons <- hg38_transcript_to_plot %>% 
      arrange(end) %>%
      dplyr::filter(type == "exon") %>%
      dplyr::filter(start %in% ref_intron$end | end %in% ref_intron$start)
    
    
    #############################
    ## PLOT
    #############################
    
    if ( hg38_transcript_to_plot %>% nrow() > 0 ) {
      
      missplicing_plot <- exons_to_plot %>%
        ggplot(aes(
          xstart = start,
          xend = end,
          y = hg38_transcript_to_plot$transcript_id %>% unique()
        )) +
        geom_range(
          data = cds_to_plot,
          fill = "#999999"
        ) +
        geom_range(
          data = utr_to_plot,
          fill = "#333333",
          height = 0.25
        ) +
        ggtranscript::geom_intron(
          data = ggtranscript::to_intron(exons_to_plot, "transcript_id"),
          aes(strand = hg38_transcript_to_plot$strand %>% unique)
        ) + 
        ggtranscript::geom_junction(
          data = df_intron,
          #aes(colour = novel_type),
          ncp = 100, 
          junction.y.max = 0.5 
        ) +
        #scale_colour_manual(breaks = c("novel_acceptor", "novel_donor"),
        #                    values = c("#C82803FF", "#23C3E4FF")) +
        theme_light() +
        ggforce::facet_zoom(xlim = c(min(ref_exons$start),
                                     max(ref_exons$end)),
                            zoom.data = zoom)  + 
        # geom_junction_label_repel(
        #   data = df_intron %>% dplyr::mutate( zoom = TRUE ) ,
        #   aes(label = paste0("seen in ", n_individuals, " out of ", db_metadata %>% nrow, " samples")),
        #   junction.y.max = 0.5
        # ) + 
        theme(axis.text.y = element_text(angle = 90, 
                                         hjust = 0.5,
                                         size = "11"),
              axis.text = element_text(size = "11"),
              axis.title = element_text(size = "11"),
              legend.position = "top",
              legend.text = element_text(size = "11"),
              legend.title = element_text(size = "11")) +
        xlab(paste0("Genomic position (chr", hg38_transcript_to_plot$seqnames %>% unique ,")")) + 
        ylab("hg38 Transcript") +
        guides(color = guide_legend(title = "Junction type: ")) +
        labs(title = paste0("Coordinates: '", df_intron$ID %>% unique, "'"))
      
      missplicing_plot %>%
        return()
      
      
    } else {
      ggplot() +
        theme_void(base_size = 14) +
        geom_text(aes(0,0,label='Sorry, it has not been possible to plot the selected junction. Please, try again.')) %>%
        return()
    }
  }
  
}

setup_ui_details_section <- function(query_results) {
  
  v <- list()
  i <- 1
  common_jxn_info <- query_results %>% distinct(coordinates, transcript_id, .keep_all = T) %>% drop_na(ref_junID)
  
  if ( str_detect(string = common_jxn_info$junID_type %>% unique(), pattern = "novel" ) ) {
    
    info_jxn_header <- c("Distance: ", "Reference intron: ")
    info_jxn <- c(paste0(common_jxn_info$distance %>% unique, " bp to annotated intron."),
                  paste0(paste(common_jxn_info$ref_coordinates %>% unique, collapse = " ")))
    if (common_jxn_info$junID_type %>% unique() == "novel_donor") {
      text_color = "#35B779FF"
    } else if ( common_jxn_info$junID_type %>% unique() == "novel_acceptor") {
      text_color = "#64037d"
    }
  } else {
    info_jxn_header <- c("","")
    info_jxn <- c("","")
    text_color = "blue"
  } 
  
  
  v[[i]] <- div(
    h3("Details"),
    hr(),
    p(strong("Coordinates:"), span(common_jxn_info$coordinates %>% unique, id = "span_junID")),
    p(strong("Junction type:", span(str_to_title(string = str_replace(string = common_jxn_info$junID_type %>% unique, pattern = "_", replacement = " ")), style = paste0("color:", text_color))),
      " (Ensembl v105)."),
    p(strong(info_jxn_header[2]), info_jxn[2]), ## Ref junction
    p(strong(info_jxn_header[1]), info_jxn[1]), ## Distance
    p(strong("Gene:"), paste(common_jxn_info$gene_name %>% unique, collapse = " ") ),
    p(strong("Transcripts:"))
  )
  i = i + 1
  
  for (transcript in (query_results$transcript_ENS %>% unique %>% sort())) {
    v[[i]] <- 
      span(tags$a(href = "#", 
                  onclick = "var ID_jxn= $('#span_junID')[0].innerText; 
Shiny.setInputValue('junID_tab1', ID_jxn);

var ID_database= $('#h4_databaseID')[0].innerText; 
Shiny.setInputValue('databaseName_tab1', ID_database);

Shiny.setInputValue('transcriptENS_tab1', this.text );
                   $('#modalVisualiseTranscript_tab1').modal('show');", 
                  transcript ), HTML('&nbsp;'))
    
    i = i + 1
  }
  
  return(v)
}


setup_ui_results_section <- function(query_results, v, i) {
  
  
  databases_queried <- query_results %>%
    distinct(database) %>%
    drop_na() %>%
    pull() 
  
  for ( j in 1:(length(databases_queried)) ){
    
    # j = 1
    
    print(i)
    
    local_query_results <- query_results %>% 
      filter(database == databases_queried[j]) %>% 
      group_by(project, cluster) %>%
      distinct( ref_junID, novel_junID, .keep_all = T) %>%
      ungroup()
    
    
    current_db_general_info <- database_equivalences %>% filter(name == local_query_results$database %>% unique)
    
    print(current_db_general_info)
    # 1. Get avg number of individuals and reads corresponding to the junx of interest
    
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
    
    
    # 2. Get message of the number clusters and projects were the junction has been per database
    
    if ( str_detect(string = local_query_results$database %>% unique(), pattern = "AD")) {
      text_n_clusters = paste0("Found in ", paste(local_query_results$cluster %>% unique(), collapse = " and "), " samples ")
      
    } else if ( str_detect(string = local_query_results$database %>% unique(), pattern = "Brainseq") ) {
      
      text_n_clusters = paste0("Found in ", paste(local_query_results$project %>% unique() %>% sort, collapse = " and "), " samples, in ", 
                               paste(local_query_results$cluster %>% unique(), collapse = " and "), " clusters")

     
      if ( any(database_metadata[[ current_db_general_info$index ]] %>% names == "sid") ) {
        database_metadata[[ current_db_general_info$index ]] <- database_metadata[[ current_db_general_info$index ]] %>%
          dplyr::select(-sample_id) %>%
          dplyr::rename(sample_id = sid)
      }
      
    } else if ( str_detect(string = local_query_results$database %>% unique(), pattern = "ENCODE|Age")) {
      text_n_clusters = paste0("Found in ", local_query_results %>% dplyr::count(project) %>% nrow(), " independent ", 
                               local_query_results$cluster_type %>% unique, ", in ", 
                               paste(local_query_results$cluster %>% unique() %>% sort, collapse = " and "), " clusters")

    } else {
      text_n_clusters = paste0("Found in ", local_query_results %>% dplyr::count(project) %>% nrow(), 
                               " independent ", local_query_results$cluster_type %>% unique, ", in ", 
                               local_query_results %>% dplyr::count(cluster) %>% nrow(), " tissues")
      
  
      if (any(database_metadata[[ current_db_general_info$index ]] %>% names == "external_id")) {
      database_metadata[[ current_db_general_info$index ]] <- database_metadata[[ current_db_general_info$index ]] %>%
        dplyr::select(-sample_id) %>%
        dplyr::rename(sample_id = external_id)
      }
    }
    
    
   
    per_individuals_message = paste0(n_individuals, " samples (", 
                                     round((n_individuals * 100)/
                                             (database_metadata[[ current_db_general_info$index ]] %>% 
                                                group_by(SRA_project) %>% 
                                                distinct(sample_id) %>% 
                                                nrow()), digits=2),
                                     "% of ", 
                                     (database_metadata[[ current_db_general_info$index ]] %>% 
                                        group_by(SRA_project) %>% 
                                        distinct(sample_id) %>% 
                                        nrow())," total samples studied)")
    
    ## Tidy database results prior returning to UI
    
    
    ## 1. Relocate junctions
    if (local_query_results$junID_type %>% unique == "annotated intron") {
      ## If it is a NEVER-MISSPLICED annotated intron, we remove data corresponding to the novel junction
      local_query_results <- local_query_results %>%
        dplyr::relocate(dplyr::any_of(x = c("ref_coordinates", "novel_coordinates")))
      
      if (local_query_results$ref_type %>% unique == "never") {
        
        local_query_results <- local_query_results %>%
          dplyr::select(-c(novel_coordinates, 
                           novel_n_individuals,
                           novel_sum_counts))
      }
      
    } else {
      
      local_query_results <- local_query_results %>%
        dplyr::relocate("novel_coordinates","novel_n_individuals", "novel_sum_counts", "novel_type", "ref_coordinates")
    }
    
    local_query_results <- local_query_results %>%
      dplyr::relocate(gene_tpm, .before = gene_name) %>%
      dplyr::relocate("cluster", "project", "database", .after = last_col()) 
    
    
    ## 2. Rename junctions
    rename_column_lookup <- c(intron_coordinates = "ref_coordinates",
                              novel_jxn_coordinates = "novel_coordinates",
                              intron_n_individuals = "ref_n_individuals",
                              intron_sum_reads = "ref_sum_counts",
                              intron_misspliced_at = "ref_type",
                              intron_donor_sequence = "ref_donor_sequence",
                              intron_acceptor_sequence = "ref_acceptor_sequence",
                              intron_MES_5ss = "ref_mes5ss",
                              intron_MES_3ss = "ref_mes3ss",
                              MSR_Donor = "MSR_D",
                              MSR_Acceptor = "MSR_A",
                              novel_jxn_n_individuals = "novel_n_individuals",
                              novel_jxn_sum_reads = "novel_sum_counts",
                              novel_jxn_type = "novel_type",
                              novel_jxn_donor_sequence = "novel_donor_sequence",
                              novel_jxn_acceptor_sequence = "novel_acceptor_sequence",
                              novel_jxn_MES_5ss = "novel_mes5ss",
                              novel_jxn_MES_3ss = "novel_mes3ss",
                              sample_cluster = "cluster")

    local_query_results <- local_query_results %>%
      dplyr::select(-c(ref_junID, 
                       dplyr::any_of(x = c("novel_junID.x", "novel_junID.y", "transcript_id", "transcript_id_list", "novel_junID",
                                           "seqnames", "start", "end", "strand", "distance", "junID", "coordinates", "cluster_index", "misspliced",
                                           "u2_intron", "protein_coding", "cluster_type", "junID_type")))) %>%
      dplyr::rename(dplyr::any_of(rename_column_lookup))
    
    
    
    
    
    v[[i]] <- div(
      h4(id="h4_databaseID", local_query_results$database %>% unique()),
      p(text_n_clusters,  style = "color:blue"),
      #p(tags$li(paste(local_query_results$project, "-", local_query_results$cluster %>% unique, collapse = ", "))),
      p(strong("Found in: "), per_individuals_message),
      p(strong("Avg. reads: "), n_reads),
      DT::datatable( data =  local_query_results %>% mutate(across(where(is.numeric), round, 5)),
                     extensions = "Buttons", 
                     options = list(
                       pageLength = 5,
                       scrollX=TRUE,
                       scrollCollapse=TRUE,
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel')
                       #initComplete = I('function(setting, json) { alert("done"); }')
                     )
      ),
      hr()
    )
    i <- i + 1
  }
  
  return(v)
}

#######################################
## CALLS AND SETTING VARIABLES
#######################################

## Needed for ggtranscript
hg38_transcripts <- readRDS( file = "./dependencies/hg38_transcripts.105.rds")

database_metadata <- get_database_metadata(database_equivalences$sqlite_file)

# gtex_table_list <- database_metadata[[ database_equivalences %>% filter["splicing_1read.sqlite"] ]] %>% distinct(cluster) %>% pull() %>% sort()
# encode_table_list <- database_metadata[[ clusters_index["ENCODE_SR_1read.sqlite"] ]] %>% distinct(SRA_project) %>% pull() %>% sort()

