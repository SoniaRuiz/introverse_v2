library(tidyverse)
library(DBI)
library(ggplot2)
library(ggtranscript)

#setwd()



database_folder <- base::file.path("database")

database_equivalences <- data.frame(index = c(1,2,3,4,5,6,7),
                                    name = c("TCGA","GTEx v8", "GTEx v8 Age", "ENCODE shRNA", "Alzheimer's Disease/Control","Parkinson's Disease/Control", "LRRK2 G2019S mutation" ),
                                    sqlite_file = c("TCGA_1read_subsampleFALSE.sqlite","splicing_1read.sqlite","splicing_1read_age.sqlite","ENCODE_SR_1read.sqlite","SRP100948_1read_subsampleFALSE.sqlite", "SRP058181_1read_subsampleFALSE.sqlite", "SRP151040_1read_subsampleFALSE.sqlite" ),
                                    project_types = c("cancer types","body sites", "body sites","RBP knockdown experiments","Alzheimer's Disease/Control", "Parkinson's Disease/Control","LRRK2 G2019S mutation" ))
                                                         


#######################################
## FUNCTIONS
#######################################


query_database <- function(chr, start, end, strand, 
                           databases_name, table_name = NULL) {
  

  # chr = 1
  # start = 108896431
  # end = 108896864
  # strand = "+"
  
  # chr = 10
  # start = 87862643
  # end = 87864453
  # strand = "+"
  
  # chr = 1
  # start=155235312
  # end= 155235462
  # strand="-"

  # chr = 1
  # start = 155235308
  # end = 155235680
  # strand="-"
  table_name = NULL
   
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
  
  # databases_name = "All"
  # databases_name = "TCGA"
  # databases_name = "GTEx v8"
  # databases_name = "ENCODE shRNA"
  # databases_name = "Alzheimer's Disease/Control"
  # databases_name = "Brainseq Squizo/Control"
  # databases_name = "Parkinson's Disease/Control"
  # databases_name = "LRRK2 G2019S mutation"
  
  # print(databases_name)
  
  if ( any(databases_name == "All") ){
    databases_sqlite_list <- database_equivalences$sqlite_file
  } else {
    databases_sqlite_list <- database_equivalences %>% filter(name == databases_name) %>% pull(sqlite_file)
  }
  
  query_results <- map_df(databases_sqlite_list, function(database_sqlite) {
   
    # database_sqlite <- databases_sqlite_list[1]
    # database_sqlite <- databases_sqlite_list[2]
    print(database_sqlite)
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
      
      # main_path<-"/home/jbrenton//Wood_ASAP_post_align_files/Final_dataset/Leafcutter"
      # 
      # pathogenic_jxn <- read_csv(file = file.path(main_path, "/Data/tdp_cryptics_junct_list.csv"))
      # pathogenic_jxn$chr<-gsub(x = pathogenic_jxn$paste_into_igv_junction, pattern = "(chr.+):.*", replacement = "\\1")
      # pathogenic_jxn$start<-gsub(x = pathogenic_jxn$paste_into_igv_junction, pattern = "chr.+:([0-9]*)-.*", replacement = "\\1")
      # pathogenic_jxn$end<-gsub(x = pathogenic_jxn$paste_into_igv_junction, pattern = "chr.+:[0-9]*-([0-9]*)", replacement = "\\1")
      # pathogenic_jxn_gr <- pathogenic_jxn %>%
      #   filter(junc_cat %in% c("novel_donor", "novel_acceptor")) %>%
      #   mutate(start = (start %>% as.integer()) +1,
      #          end = (end %>% as.integer())-1)%>%
      #   mutate(paste_into_igv_junction = paste0(chr,":", start, "-", end, ":", strand)) %>%
      #   distinct(paste_into_igv_junction, .keep_all =T) %>%
      #  filter(n_significant_data ==14 )
      
     # map_df(pathogenic_jxn_gr$paste_into_igv_junction, function(jxn_of_interest_coordinates) {
        
        print(jxn_of_interest_coordinates)
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
          query <- paste0("SELECT ref_junID, novel_junID, novel_coordinates, novel_type FROM 'novel' 
                        WHERE ref_junID='", master_information$ref_junID, "'")
          novel_coordinates <- DBI::dbGetQuery(con, query) %>% as_tibble()
          
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
          print("Not results found!")
          return(NULL)
        } 
        
      } else {
        print("Not results found!")
        return(NULL)
        
      }
      
    } else {
      print("Not results found!")
      return(NULL)
    }  
    
    
  })
  
  #DBI::dbDisconnect(con)
  query_results %>%
    as_tibble() %>%
    return()

  
  
}


setup_UI_details_section <- function(query_results) {
  
  
  v <- list()
  i <- 1
  common_jxn_info <- query_results %>% distinct(coordinates, transcript_id, .keep_all = T) %>% drop_na(ref_junID)
  
  
  # Set UI data
  
  distance_and_global_info <- set_UI_distance_information(jxn_info = common_jxn_info)
  
  transcript_span_tag <- set_UI_transcript_information(query_results)

  encori_span_tag <- set_UI_clip_information(query_results)
  
  clinvar_span_tag <- set_UI_clinvar_information(jxn_info = common_jxn_info)
  
  
  ## Return UI data
  
  v[[i]] <- div(
    
    ## Set hidden info
    p(id = "phastCons17_5ss", style = "display:none;", common_jxn_info$mean_phastCons17way5ss_100 %>% unique),
    p(id = "phastCons17_3ss", style = "display:none;", common_jxn_info$mean_phastCons17way3ss_100 %>% unique),
    p(id = "CDTS_5ss", style = "display:none;", common_jxn_info$mean_CDTS5ss_100 %>% unique),
    p(id = "CDTS_3ss", style = "display:none;", common_jxn_info$mean_CDTS3ss_100 %>% unique),
    p(id = "MES_5ss", style = "display:none;", common_jxn_info$mes5ss %>% unique),
    p(id = "MES_3ss", style = "display:none;", common_jxn_info$mes3ss %>% unique),
    p(id = "donor_sequence", style = "display:none;", common_jxn_info$donor_sequence %>% unique),
    p(id = "acceptor_sequence", style = "display:none;", common_jxn_info$acceptor_sequence %>% unique),
    
    ## Set visible info
    h3("Junction Details"),
    hr(),
    p(strong("Coordinates:"), span(common_jxn_info$coordinates %>% unique, id = "span_junID")),
    p(strong("Length (in bp):"), span(common_jxn_info$length %>% unique)),
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
    
    local_query_results <- query_results %>% 
      filter(database == databases_queried[j]) %>% 
      group_by(project, cluster) %>%
      distinct( ref_junID, novel_junID, .keep_all = T) %>%
      ungroup()
    
    
    current_db_general_info <- database_equivalences %>% filter(name == local_query_results$database %>% unique)
    
    print(current_db_general_info)
    
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
    
    # if ( str_detect(string = local_query_results$database %>% unique(), pattern = "Alzhe") ) {
    #   text_n_clusters = paste0("Found in ", paste(local_query_results$cluster %>% unique(), collapse = " and "), " samples.")
    #   
    # } else if ( str_detect(string = local_query_results$database %>% unique(), pattern = "Brainseq") ) {
    #   
    #   text_n_clusters = paste0("Found in ", paste(local_query_results$project %>% unique() %>% sort, collapse = " and "), " samples, in ", 
    #                            paste(local_query_results$cluster %>% unique(), collapse = " and "), " clusters.")
    # 
    #  
    #   if ( any(database_metadata[[ current_db_general_info$index ]] %>% names == "sid") ) {
    #     database_metadata[[ current_db_general_info$index ]] <- database_metadata[[ current_db_general_info$index ]] %>%
    #       dplyr::select(-sample_id) %>%
    #       dplyr::rename(sample_id = sid)
    #   }
    #   
    # } else if ( str_detect(string = local_query_results$database %>% unique(), pattern = "ENCODE") ) {
    #   text_n_clusters = paste0("Found in ", local_query_results %>% dplyr::count(project) %>% nrow(), " independent ", 
    #                            local_query_results$cluster_type %>% unique, ", in ", 
    #                            paste(local_query_results$cluster %>% unique() %>% sort, collapse = " and "), " samples.")
    # 
    # } else if ( str_detect(string = local_query_results$database %>% unique(), pattern = "Age") ) {
    #   text_n_clusters = paste0("Found in ", local_query_results %>% dplyr::count(project) %>% nrow(), " independent ", 
    #                            local_query_results$cluster_type %>% unique, ", in '", 
    #                            paste(local_query_results$cluster %>% unique() %>% sort, collapse = "' and '"), "' years-old samples.")
    #   
    # } else {
    #   text_n_clusters = paste0("Found in ", local_query_results %>% dplyr::count(project) %>% nrow(), 
    #                            " independent ", local_query_results$cluster_type %>% unique, ", in '", 
    #                            paste(local_query_results$cluster %>% unique, collapse = "', '"),"' samples.")
    #   
    # 
    #   if ( any(database_metadata[[ current_db_general_info$index ]] %>% names == "external_id")) {
    #     
    #     if ( any(database_metadata[[ current_db_general_info$index ]] %>% names == "sample_id") ) {
    #       database_metadata[[ current_db_general_info$index ]] <- database_metadata[[ current_db_general_info$index ]] %>%
    #         dplyr::select(-sample_id) %>%
    #         dplyr::rename(sample_id = external_id)
    #     } else{
    #       database_metadata[[ current_db_general_info$index ]] <- database_metadata[[ current_db_general_info$index ]] %>%
    #         dplyr::rename(sample_id = external_id)
    #     }
    #     
    #   }
    # }
    
    
   
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
      if ( any(local_query_results$ref_type %>% unique == "never") ) {
        
        local_query_results <- local_query_results %>%
          dplyr::select(-c(novel_coordinates, 
                           novel_n_individuals,
                           novel_sum_counts,
                           novel_type))
        
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
                                           "donor_sequence", "acceptor_sequence", "length","database","gene_name",
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
                                                    #p(text_n_clusters,  style = "color:blue"),
                                                    p(strong("Found in: "), span(per_individuals_message,  style = "color:blue")),
                                                    p(strong("Avg. reads: "), span(n_reads, style = "color:blue")),
                                                    DT::datatable( data =  final_query_results,
                                                                   rownames = F,
                                                                   extensions = c('Buttons', 'Responsive'),
                                                                   options = list(
                                                                     lengthMenu = c(5, 10, 20, 40, 80,100),
                                                                     pageLength = 5,
                                                                     scrollX=TRUE,
                                                                     scrollY=TRUE,
                                                                     scrollCollapse=TRUE,
                                                                     dom = 'Bfrtip',
                                                                     buttons = c('copy', 'csv', 'excel')
                                                                     #initComplete = I('function(setting, json) { alert("done"); }')
                                                                     ))
                                   ),
                                 hr()   
    )
  }
  #result_tabs

  return(result_tabs)
}



#######################################
## SET INITIAL VARIABLES
#######################################

## Needed for ggtranscript
# hg_111 <- rtracklayer::import(con = "/mnt/PROJECTS/splicing-accuracy-manuscript/dependencies/Homo_sapiens.GRCh38.111.chr.gtf")
# hg_111 %>%
#   as_tibble %>%
#   mutate(seqnames = paste0("chr",seqnames)) %>%
#   filter(type %in% c("transcript","exon","CDS","five_prime_utr","three_prime_utr")) %>%
#   dplyr::select(seqnames, start, end, width, type, strand, transcript_id, transcript_biotype, gene_name) %>%
#   saveRDS( file = "./dependencies/hg38_transcripts.111.rds")
# hg38_transcripts <- readRDS( file = "./dependencies/hg38_transcripts.111.rds")


## Needed for binding sites info
# encori_data <- read.delim( file = "./dependencies/ENCORI_hg38_RBPTarget_all_all_all.txt") %>%
#   dplyr::select(seqnames = chromosome, start = narrowStart, end = narrowEnd, strand = strand, RBP, geneName) %>%
#   drop_na() %>%
#   GenomicRanges::GRanges()
# MANE <- readRDS(file = "./dependencies/MANE_genes_CLIP_sites.rds")

database_metadata <- get_database_metadata(database_equivalences$sqlite_file)

# gtex_table_list <- database_metadata[[ database_equivalences %>% filter["splicing_1read.sqlite"] ]] %>% distinct(cluster) %>% pull() %>% sort()
# encode_table_list <- database_metadata[[ clusters_index["ENCODE_SR_1read.sqlite"] ]] %>% distinct(SRA_project) %>% pull() %>% sort()

