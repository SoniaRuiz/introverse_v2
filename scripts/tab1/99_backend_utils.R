## source("/home/sruiz/POST_DOC/introverse_v2/scripts/99_backend_utils.R")
library(tidyverse)

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



get_database_genes <- function(database_list) {
  
  purrr::map(database_list, function(database_name) {
    
    print(database_name)
    # database_name <- database_list[2]
    database_path <- file.path(database_folder, database_name)
    
    con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
    
    query <- paste0("SELECT DISTINCT gene_name FROM 'gene'")
    gene_name <- (DBI::dbGetQuery(con, query))$gene_name
    
    DBI::dbDisconnect(conn = con)
    
    gene_name %>%
      return()
    
  })
  
}



create_CLIP_database_tables <- function() {
  
  # database_genes <- get_database_genes(database_equivalences$sqlite_file) %>% unlist %>% unname %>% unique %>% sort()
  
  # write.table(x = database_genes, 
  #             file = "./dependencies/gene_names.txt", row.names = F, col.names = F, quote = F)
  
  
  database_genes <- readr::read_delim(file = file.path(here::here(), "dependencies/gene_names.txt"), delim = "\t", col_names = F,show_col_types = FALSE )[-1,]$X1 
  database_path <- file.path(here::here(), "database/clip_data.sqlite")
  
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  for (gene_name in database_genes) {
    
    #print(gene_name)
    
    # gene_name <- database_genes[2]
   
    if ( DBI::dbExistsTable(con, gene_name) ) {
      
      logger::log_info(paste0("'", gene_name, "' table already exists!"))
      
    } else {
      
      gene_clip_file_name <- file.path(here::here(), "dependencies/ENCORI", paste0("ENCORI_hg38_RBP_all_Gene_", gene_name, ".txt"))
      
      if ( file.exists(gene_clip_file_name) ) {
      
        df_clip_gene <- readr::read_delim(file = gene_clip_file_name, 
                                          delim = "\t", col_names = T, skip = 3, show_col_types = F)
        
        if ( df_clip_gene %>% nrow() > 0 ) {
          
          
          DBI::dbWriteTable(conn = con,
                            name = gene_name,
                            value = df_clip_gene %>%
                              dplyr::rename(seqnames = chromosome, start = narrowStart, end = narrowEnd) %>% 
                              dplyr::select(-c(geneID,clusterID,clusterNum,totalClipSiteNum,clipExpNum ,"broadStart", "broadEnd",`HepG2(shRNA)`, `K562(shRNA)`, `HepG2(CRISPR)`, `K562(CRISPR)`)),
                            overwrite = T)
          
          
          
          logger::log_info(paste0("Table: '", gene_name, "' created!"))  
          
        } else {
          logger::log_info(paste0(gene_name, " does not have CLIP data!"))
        }
        
      } else {
        logger::log_info(paste0(gene_name, " does not have files downloaded from ENCORI!"))
      }
    
    }
    
    
  }
  DBI::dbDisconnect(conn = con)
}



create_hg38_transcript_table <- function() {
  
  MANE <- readRDS(file = "./dependencies/MANE_genes_CLIP_sites.rds")
  hg38_transcripts <- readRDS( file = "./dependencies/hg38_transcripts.111.rds")
  
  
  
  
  database_path <- file.path(here::here(), "database/hg38_transcripts.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  
  hg38_transcripts_tidy <- hg38_transcripts %>%
    left_join(y = MANE %>% 
                dplyr::select(transcript_id, tag) %>% 
                mutate(transcript_id = stringr::str_sub(start = 1, end = 15, string = transcript_id)) %>%
                distinct(transcript_id, .keep_all = T),
              by = "transcript_id") %>%
    mutate(type = type %>% as.character()) %>%
    #rowwise() %>%
    mutate(type = ifelse(str_detect(type, pattern = "utr"), "UTR", type)) %>%
    ungroup()
  
  hg38_transcripts_tidy$tag %>% unique
  
  DBI::dbWriteTable(conn = con,
                    name = "hg38_transcripts",
                    value = hg38_transcripts_tidy,
                    overwrite = T)
  
  DBI::dbDisconnect(conn = con)
  
  # hg_111 <- rtracklayer::import(con = "/mnt/PROJECTS/splicing-accuracy-manuscript/dependencies/Homo_sapiens.GRCh38.111.chr.gtf")
  # hg_111 %>%
  #   as_tibble %>%
  #   mutate(seqnames = paste0("chr",seqnames)) %>%
  #   filter(type %in% c("transcript","exon","CDS","five_prime_utr","three_prime_utr")) %>%
  #   dplyr::select(seqnames, start, end, width, type, strand, transcript_id, transcript_biotype, gene_name) %>%
  #   saveRDS( file = "./dependencies/hg38_transcripts.111.rds")
  
}



get_genomic_coordinates <- function(coordinates) {
  
  if (str_detect(coordinates, pattern = ",")) {
    coordinates <- str_split(string = coordinates, pattern = ",") %>% unlist
  }
  
  map_df(coordinates, function(coordinate) {
    
    # coordinate <- df_gene_splicing$novel_coordinates[1]
    
    chr_junc <- coordinate %>%
      str_sub(start = 1,
              end = str_locate_all(string = coordinate, pattern = ":")[[1]][1,2]-1)
    
    ## Clinvar mutation
    if (str_locate_all(string = coordinate, pattern = "-")[[1]][,1] %>% length() == 0) {
      
      start_junc <- coordinate %>%
        str_sub(start = str_locate_all(string = coordinate, pattern = ":")[[1]][1,2]+1,
                end = coordinate %>% stringr::str_count() )
      end_junc <- start_junc
      
      data.frame(ID = coordinate,
                 seqnames = chr_junc,
                 start = start_junc %>% as.integer(),
                 end = end_junc %>% as.integer()) %>%
        return()
      
    } else {
      start_junc <- coordinate %>%
        str_sub(start = str_locate_all(string = coordinate, pattern = ":")[[1]][1,2]+1,
                end = str_locate_all(string = coordinate, pattern = "-")[[1]][1,2]-1)
      end_junc <- coordinate %>%
        str_sub(start = str_locate_all(string = coordinate, pattern = "-")[[1]][1,2]+1,
                end = str_locate_all(string = coordinate, pattern = ":")[[1]][2,2]-1)
      strand_junc <- coordinate %>%
        str_sub(start = str_locate_all(string = coordinate, pattern = ":")[[1]][2,2]+1,
                end = coordinate %>% stringr::str_count())
      
      data.frame(ID = coordinate,
                 seqnames = chr_junc,
                 start = start_junc %>% as.integer(),
                 end = end_junc %>% as.integer(),
                 strand = strand_junc) %>%
        return()
    }
    
    
  })
}

###############################################
## CALLS 
###############################################
# 
# create_CLIP_database_tables()
