## source("/home/sruiz/POST_DOC/introverse_v2/scripts/99_backend_utils.R")

database_connect <- function(database_sqlite_list, project = NULL) {
  
  conn_list <- list()
  for (database_sqlite in database_sqlite_list) {
    
    # database_sqlite <- database_sqlite_list[2]
    database_path <- file.path(database_folder, database_sqlite)
    
    con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
    
    conn_list[database_sqlite] = con
    
  }
  return(conn_list)
  
}



get_database_metadata <- function(database_sqlite_list, project = NULL) {
  
  purrr::map(database_sqlite_list, function(database_sqlite) {
    
    # database_sqlite <- database_sqlite_list[2]
    # database_sqlite <- "TCGA_1read_subsampleFALSE.sqlite"
    # project <- "BLCA"
    
    # database_path <- file.path(database_folder, database_sqlite)
    # con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
    con <- conn_list[database_sqlite][[1]]
    
    query <- paste0("SELECT * FROM 'metadata'")
    metadata <- DBI::dbGetQuery(con, query)
    #DBI::dbDisconnect(conn = con)
    
    if ( !is.null(project) ) {
      metadata <- metadata %>%
        filter(SRA_project == project)
    }
    
    metadata %>%
      return()
    
  })
  
}



get_database_genes <- function(conn_list) {
  
  purrr::map(conn_list, function(con) {
    
    query <- paste0("SELECT gene_id, gene_name FROM 'gene'")
    (DBI::dbGetQuery(con, query)) %>%
      as_tibble() %>%
      drop_na() %>%
      arrange(gene_name) %>%
      return()
    
  })
  
}



get_database_transcripts <- function(con, gene_name) {
  
 # gene_name = "PSEN1"
    
  query <- paste0("SELECT gene.gene_name, transcript.transcript_id, transcript.MANE 
                  FROM 'transcript' 
                  INNER JOIN 'gene' ON transcript.gene_id = gene.id
                  WHERE gene.gene_name = '", gene_name, "'")
  transcripts_to_return <- (DBI::dbGetQuery(con, query)) %>%
    as_tibble() %>%
    drop_na() %>%
    arrange(desc(MANE), transcript_id)
  
  
  database_path <- file.path(here::here(), "database/hg38_transcripts.sqlite")
  cont <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbListTables(cont)
  transcripts_w_biotype <- DBI::dbGetQuery(conn = cont, statement = paste0("SELECT DISTINCT transcript_id, transcript_biotype FROM'hg38_transcripts' 
                                                  WHERE transcript_id IN ('",paste(transcripts_to_return$transcript_id, collapse="','"),"')"))

  transcripts_to_return %>%
    inner_join(transcripts_w_biotype,
               by = "transcript_id") %>%
    as_tibble() %>%
    drop_na() %>%
    arrange(desc(MANE), transcript_id) %>%
    mutate(transcript_label = ifelse(MANE, paste0(transcript_id, " (MANE)"), 
                                     paste0(transcript_id, " (",str_replace_all(string = transcript_biotype, pattern = "_", replacement = " "),")"))) %>%
    return()
    

  
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
  #DBI::dbDisconnect(conn = con)
}



create_hg38_database_table <- function() {
  
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
  
  #DBI::dbDisconnect(conn = con)
  
  # hg_111 <- rtracklayer::import(con = "/mnt/PROJECTS/splicing-accuracy-manuscript/dependencies/Homo_sapiens.GRCh38.111.chr.gtf")
  # hg_111 %>%
  #   as_tibble %>%
  #   mutate(seqnames = paste0("chr",seqnames)) %>%
  #   filter(type %in% c("transcript","exon","CDS","five_prime_utr","three_prime_utr")) %>%
  #   dplyr::select(seqnames, start, end, width, type, strand, transcript_id, transcript_biotype, gene_name) %>%
  #   saveRDS( file = "./dependencies/hg38_transcripts.111.rds")
  
}



create_BigWig_URL_list <- function() {
  
  
  packageVersion("recount3")
  library("recount3")
  library("tidyverse")
  
  
  human_projects <- recount3::available_projects(organism = "human")
  
  ## GET RECOUNT3 BIGWIG FILES
  bigWig_URLs_recount3 <- map_df(c("tcga","gtex","sra"), function(source) {
    
    recount_projects <- human_projects[human_projects$file_source == source, ] %>%
      as_tibble()
    
    # source = "sra"
    if (source == "sra") {
      recount_projects <- recount_projects %>%
        filter(project == "SRP100948")
    }
    
    
    map_df(recount_projects$project, function(project_id) {
      
      # project_id <- recount_projects$project[1]
      message(project_id, "...")
      
      if (source == "gtex") {
        base_URL = "~/PROJECTS/splicing-accuracy-manuscript/results/splicing_1read/111/"
      } else if (source == "sra") {
        base_URL = "~/PROJECTS/recount3-database-project/results/SRP100948_1read_subsampleFALSE/111/"
      } else {
        base_URL = "~/PROJECTS/recount3-database-project/results/TCGA_1read_subsampleFALSE/111/"
      }
      
      # https://recount-opendata.s3.amazonaws.com/recount3/release
      # https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/gtex/base_sums/EL/BLOOD_VESSEL/DY/gtex.base_sums.BLOOD_VESSEL_GTEX-14BMV-1226-SM-5TDDY.1.ALL.bw
      metadata_path = file.path(base_URL, project_id, "base_data", paste0(project_id, "_samples_raw_metadata.rds"))
      
      if ( file.exists(metadata_path) ) {
        
        metadata <- readRDS(file = metadata_path) %>% as_tibble() %>%
          rowwise() %>%
          mutate(BigWigURL = str_replace_all(string = BigWigURL, pattern = "http://duffel.rail.bio/recount3/", replacement = "https://recount-opendata.s3.amazonaws.com/recount3/release/")) %>%
          ungroup()
        
        if (source == "gtex") {
          
          data.frame(project = project_id,
                     cluster = metadata$gtex.smtsd,
                     BigWigURL = metadata$BigWigURL) %>%
            as_tibble()%>%
            return()
          
        } else  if (source == "sra") {
          
          metadata_tidy <- metadata %>%
            as_tibble() %>%
            dplyr::select(external_id, 
                          sra.experiment_title, 
                          sra.sample_attributes, 
                          all_mapped_reads = recount_qc.star.all_mapped_reads,
                          BigWigURL) %>%
            mutate(rn = row_number()) %>%
            separate_rows(sra.sample_attributes, sep = "\\|\\s*") %>%
            separate(sra.sample_attributes, into = c('col1', 'col2'), sep = ";;") %>% 
            pivot_wider(names_from = col1, values_from = col2) %>% 
            dplyr::select(-rn) %>%
            as.data.frame() %>%
            dplyr::select(-any_of("sample_id")) %>%
            dplyr::rename(sample_id = external_id) %>%
            mutate(diagnosis = str_remove_all(diagnosis,pattern = "'"))
          
          data.frame(project = project_id,
                     cluster = metadata_tidy$diagnosis %>% as.character(),
                     BigWigURL = metadata_tidy$BigWigURL) %>%
            return()
        } else {
          
          metadata <- metadata %>%
            drop_na(tcga.cgc_sample_sample_type)
          
          if ( metadata$tcga.cgc_sample_sample_type %>% unique %>% length() > 1) {
            
            matched_samples <- metadata %>%
              dplyr::count(tcga.gdc_cases.case_id ) %>%
              filter(n>=2) %>%
              pull(tcga.gdc_cases.case_id)
            
            
            metadata <- metadata %>%
              filter(tcga.gdc_cases.case_id %in% matched_samples)
            
            data.frame(project = project_id,
                       cluster = metadata$tcga.cgc_sample_sample_type %>% as.character(),
                       BigWigURL = metadata$BigWigURL) %>%
              return()
            
          } else {
            
            data.frame(project = project_id,
                       cluster = metadata$tcga.cgc_sample_sample_type %>% as.character(),
                       BigWigURL = metadata$BigWigURL) %>%
              return()
          }
          
          
        }
      } else {
        return(NULL)
      }
      
    })
  })
  
  
  ## GET ENCODE BIGWIG FILES
  ENCODE_metadata_path = file.path("~/PROJECTS/ENCODE_Metadata_Extraction/results/metadata_shRNA_bigWig_samples.tsv")
  ENCODE_metadata <- read.delim(file = ENCODE_metadata_path) %>% as_tibble()
  bigWig_URLs_ENCODE <- map_df(ENCODE_metadata$target_gene %>% unique, function(RBP) {
    
    message(RBP, "...")
    # RBP = (ENCODE_metadata$target_gene %>% unique)[1]
    data.frame(project = RBP,
               cluster = ENCODE_metadata %>%
                 filter(target_gene == RBP) %>%
                 distinct(sample_id, .keep_all=T) %>%
                 pull(experiment_type) %>%
                 str_to_title(),
               BigWigURL = ENCODE_metadata %>%
                 filter(target_gene == RBP) %>%
                 distinct(sample_id, .keep_all=T) %>%
                 mutate(URL = paste0("https://www.encodeproject.org/files/",sample_id,"/@@download/",sample_id,".bigWig")) %>%
                 dplyr::pull(URL),
               BigWig_type = ENCODE_metadata %>%
                 filter(target_gene == RBP) %>%
                 distinct(sample_id, .keep_all=T) %>%
                 mutate(output_type = str_remove_all(string = output_type, pattern = " strand signal of unique reads")) %>%
                 dplyr::pull(output_type) ) %>%
      as_tibble()%>%
      return()
  })
  
  
  ## JOIN BIGWIG DATASETS
  bigWig_URLs = plyr::rbind.fill(bigWig_URLs_recount3 %>% filter(cluster != "" ) %>% drop_na(), bigWig_URLs_ENCODE) %>% as_tibble()
  
  
  ## SAVE DATA
  saveRDS(object = bigWig_URLs, file = "~/POST_DOC/introverse_v2/dependencies/bigWig_URLs.rds")
  
  
  ## STORE BIGWIG DATA IN DATABASE
  # bigWig_URLs <- readRDS(file = file.path(here::here(), "dependencies/bigWig_URLs.rds")) %>% filter(str_detect(BigWigURL, pattern = "encode", negate = T))
  # 
  saveRDS(object = plyr::rbind.fill(bigWig_URLs, bigWig_URLs_ENCODE),
          file = "/mnt/POST_DOC/introverse_v2/dependencies/bigWig_URLs.rds")
  
  # database_bigwig_data <- readRDS(file = "dependencies/bigWig_URLs.rds") %>%
  #   dplyr::rename(bigwig = BigWigURL) 
  # database_path <- file.path(here::here(), "database/bigwig.sqlite")
  # 
  # con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  # 
  # DBI::dbWriteTable(conn = con,
  #                   name = "bigwig",
  #                   value = database_bigwig_data,
  #                   overwrite = T)
  
  #DBI::dbDisconnect(conn = con)
}


#########################################
## UTILS TO TRANSCRIPT VISUALISATION
#########################################

get_transcript_to_plot <- function(junID = NULL,
                                   geneName = NULL,
                                   transcript_id = NULL,
                                   multiple = F) {
  
  
  ## Connect to the 'hg38 transcripts' data database and retrieve data from the current gene
  database_path <- file.path(here::here(), "/database/hg38_transcripts.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  if ( !is.null(geneName) && !is.null(transcript_id) ) {
    query <- paste0("SELECT * FROM 'hg38_transcripts' WHERE gene_name='", geneName, "' AND transcript_id!='", transcript_id, "'")  
  } else if ( !is.null(geneName) ) {
    query <- paste0("SELECT * FROM 'hg38_transcripts' WHERE gene_name='", geneName, "'")  
  } else {
    query <- paste0("SELECT * FROM 'hg38_transcripts' WHERE transcript_id='", transcript_id, "'")  
  }
  db_transcripts_data <- DBI::dbGetQuery(con, query) 
  DBI::dbDisconnect(conn = con)
  
  
  if ( !is.null(junID) ) {
    junction_to_plot <- get_genomic_coordinates(junID)  
    db_transcripts_data <- db_transcripts_data %>%
      filter(transcript_id %in% (db_transcripts_data %>%
                                   filter(type == "transcript") %>%
                                   filter(start < junction_to_plot$start,
                                          end > junction_to_plot$end) %>%
                                   pull(transcript_id)))
  }
  
  
  ## We get the first protein-coding transcript
  if ( !multiple ) {
    
    db_transcript_data <- db_transcripts_data %>%
      filter(transcript_id == (db_transcripts_data %>%
                                 filter( transcript_biotype == "protein_coding") %>%
                                 arrange(transcript_id) %>%
                                 pull(transcript_id) %>%
                                 unique)[1] ) 
    
  } else {
    db_transcript_data <- db_transcripts_data
  }
  
  
  
  
  ## GET EXON COORDINATES FROM MANE TRANSCRIPT
  transcript <- db_transcript_data %>% dplyr::filter(type == "transcript") %>% as_tibble()
  exons <- db_transcript_data %>% dplyr::filter(type == "exon") %>% as_tibble()
  CDS <- db_transcript_data %>% dplyr::filter(type == "CDS") %>% as_tibble()
  UTR <- db_transcript_data %>% dplyr::filter(type == "UTR") %>% as_tibble()
  
  
  if ( multiple ) {
    
    if ( CDS %>% filter(transcript_biotype == "protein_coding") %>% nrow() > 0 &&
         UTR  %>% filter(transcript_biotype == "protein_coding") %>% nrow() > 0 ) {
      transcripts_to_plot <- intersect(CDS %>% filter(transcript_biotype == "protein_coding") %>% pull(transcript_id) %>% unique,
                                       UTR  %>% filter(transcript_biotype == "protein_coding") %>% pull(transcript_id) %>% unique)
      
    } else {
      transcripts_to_plot <- transcript$transcript_id %>% unique
    }
    
    exons <- exons %>%
      filter(transcript_id %in% transcripts_to_plot)
    CDS <- CDS %>%
      filter(transcript_id %in% transcripts_to_plot)
    UTR <- UTR %>%
      filter(transcript_id %in% transcripts_to_plot)
    
  }
  
  if ( exons %>% nrow > 0 | CDS %>% nrow > 0 | UTR %>% nrow > 0 ) {
    return(list(transcript = transcript,
                exons = exons,
                cds = CDS,
                utr = UTR))
  } else {
    
    return(NULL)
  }
} 


visualise_empty_plot <- function(text = NULL) {
  
  if (is.null(text)) {
    text = paste0("There are no protein-coding transcript structures available for the selected gene.")
  }
  
  ggplot() +
    theme_void() +
    geom_text(aes(x = 0, y = 0, label=text), size = 6) %>%
    return()
  
}


get_genomic_coordinates <- function(coordinates) {
  
  if (str_detect(coordinates, pattern = "; ")) {
    coordinates <- str_split(string = coordinates, pattern = "; ") %>% unlist
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