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



get_database_clinvar <- function(jxn_info) {
  
  ## Connect to the CLIP data database and retrieve data from the current gene
  database_path <- file.path("./database/clinvar.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  query <- paste0("SELECT * FROM 'clinvar' WHERE seqnames=",str_remove(string = jxn_info %>% distinct(seqnames) %>% pull(), pattern = "chr"))
  db_clip_data_gr <- DBI::dbGetQuery(con, query) %>%
    GenomicRanges::GRanges()
  DBI::dbDisconnect(conn = con)
  
  ## Find overlaps between clinvar mutations and mis-splicing ratios
  overlaps <- GenomicRanges::findOverlaps(query = db_clip_data_gr,
                                          subject = jxn_info %>% mutate(seqnames = str_remove(string = seqnames, pattern = "chr")) %>% GenomicRanges::GRanges(),
                                          
                                          type = "any")
  
  db_clip_data_gr[S4Vectors::queryHits(overlaps), ] %>%
    as_tibble() %>%
    dplyr::relocate(ID) %>%
    return()
}



create_indexes_database <- function() {
  
  ################################
  ## SPLICING DATABASES
  ################################
 
  
  for (con in conn_list) {
    
    ## INTRON TABLE -----------------------------------------------------------------
    
    ## verify indexes exist on 'intron' master table
    query <- paste0("SELECT * FROM 'sqlite_master' 
                  WHERE tbl_name = 'intron'
                  AND name = 'index_intron_coord'")
    
    
    if ( nrow(dbGetQuery(conn = con, query)) == 0) {
      
      message("Creating index 'index_intron_coord'...")
      query <- paste0("CREATE UNIQUE INDEX 'index_intron_coord' ON 'intron'(ref_coordinates)");
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
    }
    
    
    
    ## NOVEL TABLE -----------------------------------------------------------------
    
    
    ## verify indexes exist on 'novel' master table
    query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'novel' AND name = 'index_novel'")
    
    if ( nrow(dbGetQuery(conn = con, query)) == 0 ) {
      
      message("Creating index 'index_novel'...")
      query <- paste0("CREATE UNIQUE INDEX 'index_novel' ON 'novel'(ref_junID,novel_junID)");
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
    }
    
    
    ## verify indexes exist on 'novel' master table
    query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'novel' AND name = 'index_novel_coord'")
    
    if ( nrow(dbGetQuery(conn = con, query)) == 0 ) {
      
      message("Creating index 'index_novel_coord'...")
      query <- paste0("CREATE UNIQUE INDEX 'index_novel_coord' ON 'novel'(novel_coordinates)");
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
    }
    
    
    ## TRANSCRIPT TABLE -----------------------------------------------------------------
    
    query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'transcript' AND name = 'index_transcript_id'")
    
    if ( nrow(dbGetQuery(conn = con, query)) == 0 ) {
      
      message("Creating index 'index_transcript_id'...")
      query <- paste0("CREATE UNIQUE INDEX 'index_transcript_id' ON 'transcript'(id)");
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
    }
    
    
    query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'transcript' AND name = 'index_transcript_ensembl_id'")
    
    if ( nrow(dbGetQuery(conn = con, query)) == 0 ) {
      
      message("Creating index 'index_transcript_ensembl_id'...")
      query <- paste0("CREATE UNIQUE INDEX 'index_transcript_ensembl_id' ON 'transcript'(transcript_id)");
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
    }
    
    
    ## GENE TABLE -----------------------------------------------------------------
    
    ## verify indexes exist on 'gene' master tableg
    query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'gene' AND name = 'index_gene_id'")
    
    if ( nrow(dbGetQuery(conn = con, query)) == 0 ) {
      
      message("Creating index 'index_gene_id'...")
      query <- paste0("CREATE UNIQUE INDEX 'index_gene_id' ON 'gene'(id)");
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
    }
    
    ## verify indexes exist on 'gene' master tableg
    query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'gene' AND name = 'index_gene_ensembl_id'")
    
    if ( nrow(dbGetQuery(conn = con, query)) == 0 ) {
      
      message("Creating index 'index_gene_ensembl_id'...")
      query <- paste0("CREATE UNIQUE INDEX 'index_gene_ensembl_id' ON 'gene'(gene_id)");
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
    }
  }
  
  
  
  ################################
  ## OTHER DATABASES
  ################################
  
  
  ## bigwig DATABASE ---------------------------------------------------------------
  
  database_path <- file.path("./database/bigwig.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  query <- paste0("SELECT * FROM 'sqlite_master' 
                  WHERE tbl_name = 'bigwig'
                  AND name = 'index_bigwig'")
  
  if ( nrow(dbGetQuery(conn = con, query)) == 0) {
    
    message("Creating index 'index_bigwig'...")
    query <- paste0("CREATE UNIQUE INDEX 'index_bigwig' ON 'bigwig'(project, cluster, BigWigURL )");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
  }
  DBI::dbDisconnect(conn = con)
  
  ## clinvar DATABASE --------------------------------------------------------------
  
  database_path <- file.path("./database/clinvar.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  query <- paste0("SELECT * FROM 'sqlite_master' 
                  WHERE tbl_name = 'clinvar'
                  AND name = 'index_clinvar'")
  
  if ( nrow(dbGetQuery(conn = con, query)) == 0) {
    
    message("Creating index 'index_clinvar'...")
    query <- paste0("CREATE UNIQUE INDEX 'index_clinvar' ON 'clinvar'(seqnames,start,end)");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
  }
  DBI::dbDisconnect(conn = con)
  
  ## CLIP_ENCORI DATABASE ----------------------------------------------------------
  
  database_path <- file.path("./database/clip_data.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  tables <- DBI::dbListTables(conn = con)
  
  for (table in tables) {
    # table = tables[1]
    # table = "MATR3"
    
    query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = '",table,"' AND name = 'index_",table,"_coord'")
    if ( nrow(dbGetQuery(conn = con, query)) == 0 ) {
      
      table_data <- dbGetQuery(conn = con, paste0("SELECT * FROM '",table,"'") )
    
      if (table_data %>% dplyr::count(seqnames,start,end,RBP) %>% filter(n>1) %>% nrow() == 0) {

        message(table, " - creating index...")
        query <- paste0("CREATE UNIQUE INDEX 'index_",table,"_coord' ON '",table,"'(seqnames,start,end,RBP)");
        res <- DBI::dbSendQuery(conn = con, statement = query)
        DBI::dbClearResult(res)
        
      }else {
        message(table, " - index already created!")
      }
    } else {
      message(table, " - index already created!")
    }
    
  }
  
  DBI::dbDisconnect(conn = con)
  
  ## hg38 DATABASE -----------------------------------------------------------------
  
  database_path <- file.path("./database/hg38_transcripts.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  query <- paste0("SELECT * FROM 'sqlite_master' 
                  WHERE tbl_name = 'hg38_transcripts'
                  AND name = 'index_hg38_transcripts'")
  
  if ( nrow(dbGetQuery(conn = con, query)) == 0) {
    
    message("Creating index 'index_hg38_transcripts'...")
    query <- paste0("CREATE UNIQUE INDEX 'index_hg38_transcripts' ON 'hg38_transcripts'(gene_name,transcript_id,transcript_biotype)");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
  }
  
  DBI::dbDisconnect(conn = con)
}



create_clinvar_database <- function() {
  
  clinvar_data <- readRDS(file = "./dependencies/clinvar_splicing_pathogenic.rds") %>%
    dplyr::select(-c(ORIGIN, RS, AF_EXAC, CLNVI, AF_ESP, AF_TGP, CLNDISDBINCL, CLNDNINCL,
                     DBVARID, SCI, SCIDISDB, SCIDN, SCIREVSTAT, ONC, ONCDISDB, ONCDN, ONCREVSTAT, QUAL,FILTER)) %>%
    dplyr::relocate(end, .after = "start") %>%
    dplyr::relocate(ID) %>% 
    distinct(seqnames,start,end, .keep_all = T)
  

  database_path <- file.path(here::here(), "database/clinvar.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbWriteTable(conn = con,
                    name = "clinvar",
                    value = clinvar_data,
                    overwrite = T)
  query <- paste0("CREATE UNIQUE INDEX 'index_clinvar' ON 'clinvar'(seqnames,start,end)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  DBI::dbDisconnect(conn = con)
  
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
  
  if ( !file.exists(file.path(here::here(), "/dependencies/bigWig_URLs.rds")) ) {
    
    human_projects <- recount3::available_projects(organism = "human")
    
    ## GET RECOUNT3 BIGWIG FILES
    bigWig_URLs_recount3 <- map_df(c("tcga","gtex","sra"), function(source) {
      
      # source = "sra"
      
      recount_projects <- human_projects[human_projects$file_source == source, ] %>%
        as_tibble()
      
      
      if (source == "sra") {
        recount_projects <- recount_projects %>%
          filter(project %in% c("SRP058181","SRP181886"))
      }
      
      
      map_df(recount_projects$project, function(project_id) {
        
        # project_id <- recount_projects$project[1]
        message(project_id, "...")
        
        if (source == "gtex") {
          base_URL = "~/PROJECTS/splicing-accuracy-manuscript/results/splicing_1read/111/"
        } else if (source == "sra") {
          base_URL = paste0("~/PROJECTS/recount3-database-project/results/",project_id,"_1read_subsampleFALSE/111/")
          
        } else {
          base_URL = "~/PROJECTS/recount3-database-project/results/TCGA_1read_subsampleFALSE/111/"
        }
        
        # https://recount-opendata.s3.amazonaws.com/recount3/release
        # https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/gtex/base_sums/EL/BLOOD_VESSEL/DY/gtex.base_sums.BLOOD_VESSEL_GTEX-14BMV-1226-SM-5TDDY.1.ALL.bw
        
        metadata_path = file.path(base_URL, project_id, "base_data", paste0(project_id, "_samples_raw_metadata.rds"))
        
        if ( file.exists(metadata_path) ) {
          
          metadata <- readRDS(file = metadata_path) %>% as_tibble() %>%
            rowwise() %>%
            mutate(BigWigURL = str_replace_all(string = BigWigURL, 
                                               pattern = "http://duffel.rail.bio/recount3/", 
                                               replacement = "https://recount-opendata.s3.amazonaws.com/recount3/release/")) %>%
            ungroup()
          
          if (source == "gtex") {
            
            data.frame(project = project_id,
                       cluster = metadata$gtex.smtsd,
                       BigWigURL = metadata$BigWigURL) %>%
              as_tibble()%>%
              return()
            
          } else if (source == "sra") {
            
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
              dplyr::rename(sample_id = external_id) 
            
            if (project_id == "SRP181886") {
              
              metadata_tidy <- metadata_tidy %>%
                mutate(diagnosis = str_remove_all(diagnosis,pattern = "'"))
              
            } else {
              metadata_tidy <- metadata_tidy %>%
                mutate(diagnosis = ifelse( test = str_detect(sra.experiment_title, pattern="P"),
                                           yes = "Parkinsons Disease",
                                           no = "Control"))
            }
              
            
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
    bigWig_URLs_ENCODE <- map_df(c("crispr", "shRNA"), function(encode_proyect) {
      
      # encode_proyect = "crispr"
      ENCODE_metadata_path = file.path(paste0("~/PROJECTS/splicing-accuracy-manuscript/ENCODE_SR/ENCODE_Splicing_Analysis/metadata/metadata_",encode_proyect,"_samples.tsv"))
      ENCODE_metadata <- read.delim(file = ENCODE_metadata_path) %>% as_tibble()
      
      map_df(ENCODE_metadata$target_gene %>% unique, function(RBP) {
        
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
    })

    
    
    ## JOIN BIGWIG DATASETS
    bigWig_URLs = plyr::rbind.fill(bigWig_URLs_recount3 %>% filter(cluster != "" ) %>% drop_na(), bigWig_URLs_ENCODE) %>% as_tibble()
    
    
    ## SAVE DATA
    saveRDS(object = bigWig_URLs, file = "~/POST_DOC/introverse_v2/dependencies/bigWig_URLs.rds")
    
  } else {
    
    bigWig_URLs <- readRDS(file = file.path(here::here(), "/dependencies/bigWig_URLs.rds"))
  }
  
  
  ## STORE BIGWIG DATA IN DATABASE
  # bigWig_URLs <- readRDS(file = file.path(here::here(), "dependencies/bigWig_URLs.rds")) %>% filter(str_detect(BigWigURL, pattern = "encode", negate = T))
  # 
  # saveRDS(object = plyr::rbind.fill(bigWig_URLs, bigWig_URLs_ENCODE),
  #         file = "/mnt/POST_DOC/introverse_v2/dependencies/bigWig_URLs.rds")
  
  # database_bigwig_data <- readRDS(file = "dependencies/bigWig_URLs.rds") %>%
  #   dplyr::rename(bigwig = BigWigURL)
  
  database_path <- file.path(here::here(), "database/bigwig.sqlite")

  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)

  DBI::dbWriteTable(conn = con,
                    name = "bigwig",
                    value = bigWig_URLs,
                    overwrite = T)
  
  DBI::dbDisconnect(conn = con)
}


#########################################
## UTILS TO TRANSCRIPT VISUALISATION
#########################################



get_transcript_to_plot <- function(junID = NULL,
                                   geneName = NULL,
                                   transcript.id = NULL,
                                   jxn.type,
                                   multiple = F) {
  
  # junID = "chr1:155235312-155235462:-"
  # geneName = NULL
  # transcript.id = "ENST00000368373"
  # jxn.type = "Novel Donor"
  # multiple = F
  
  # junID = "chr19:4491836-4492014:+"
  # geneName = "ENST00000616600"
  # transcript.id = "ENST00000616600"
  # jxn.type = "Novel Acceptor"
  # multiple = F
  
  message("'get_transcript_to_plot': '", junID, "' '", geneName, "' '", transcript.id, "' '", jxn.type, "' '", multiple)
  
  ## Connect to the 'hg38 transcripts' data database and retrieve data from the current gene/transcript
  database_path <- file.path(here::here(), "/database/hg38_transcripts.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  if ( !is.null(geneName) && !is.null(transcript.id) ) {
    query <- paste0("SELECT * FROM 'hg38_transcripts' WHERE gene_name='", geneName, "' AND transcript_id!='", transcript.id, "'")  
  } else if ( !is.null(geneName) ) {
    query <- paste0("SELECT * FROM 'hg38_transcripts' WHERE gene_name='", geneName, "'")  
  } else {
    query <- paste0("SELECT * FROM 'hg38_transcripts' WHERE transcript_id='", transcript.id, "'")  
  }
  db_transcripts_data <- DBI::dbGetQuery(con, query) 
  DBI::dbDisconnect(conn = con)
  
  
  
  if ( !is.null(junID) ) {
    
    junction_to_plot <- get_genomic_coordinates(junID) 
    

    if ( jxn.type == "Annotated Intron" ) {
      
      db_transcripts_data <- db_transcripts_data %>%
        filter(transcript_id %in% (ggtranscript::to_intron(db_transcripts_data %>% filter(type == "exon"), "transcript_id") %>%
                                     mutate(end = end - 1,
                                            start = start + 1) %>%
                                     filter(end == junction_to_plot$end |
                                            start == junction_to_plot$start) %>%
                                     pull(transcript_id)))
    } else {
      
      ## To avoid plotting novel junctions with a novel splice site very far away from the limits of the transcript
      db_transcripts_data <- db_transcripts_data %>%
        filter(transcript_id %in% (db_transcripts_data %>%
                                     filter(type == "transcript",
                                            start < junction_to_plot$start |
                                            end > junction_to_plot$end) %>%
                                     pull(transcript_id))) %>%
        filter(transcript_id %in% (ggtranscript::to_intron(db_transcripts_data %>% filter(type == "exon"), "transcript_id") %>%
                                     mutate(end = end - 1,
                                            start = start + 1) %>%
                                     filter(end == junction_to_plot$end | 
                                              start == junction_to_plot$start) %>%
                                     pull(transcript_id)))
    }
    
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


get_exons_to_zoom <- function(jun.type, 
                              intron.to.zoom,
                              transcript.to.plot) {
  
  #message("get_exons_to_zoom: ",jun.type, " - ", intron.to.zoom, " - ", transcript.to.plot)
  
  if ( (jun.type == "Novel Donor" && (transcript.to.plot$exons$strand %>% unique) == "+") | 
       (jun.type == "Novel Acceptor" && (transcript.to.plot$exons$strand %>% unique) == "-") ) {
    
    index_first_exon <- which(abs(transcript.to.plot$exons$start - intron.to.zoom$end) == 
                                min(abs(transcript.to.plot$exons$start - intron.to.zoom$end)))
    
    index_second_exon <- which((abs(transcript.to.plot$exons$end - intron.to.zoom$start) == 
                                  min(abs(transcript.to.plot$exons$end - intron.to.zoom$start)[c(1:(index_first_exon-1))])))
    
  } else if ( (jun.type == "Novel Acceptor" && (transcript.to.plot$exons$strand %>% unique) == "+") | 
              (jun.type == "Novel Donor" && (transcript.to.plot$exons$strand %>% unique) == "-") ) {
    
    index_first_exon <- which(abs(transcript.to.plot$exons$end - intron.to.zoom$start) == 
                                min(abs(transcript.to.plot$exons$end - intron.to.zoom$start)))
    
    if (jun.type == "Novel Donor") {
      index_second_exon <- which((abs(transcript.to.plot$exons$start - intron.to.zoom$end) == 
                                    min(abs(transcript.to.plot$exons$start - intron.to.zoom$end)[-index_first_exon])))
    } else {
      index_second_exon <- which((abs(transcript.to.plot$exons$start - intron.to.zoom$end) == 
                                    min(abs(transcript.to.plot$exons$start - intron.to.zoom$end)[-c(1:(index_first_exon))])))
    }
    
    
  } else {
    
    index_first_exon <- which(abs(transcript.to.plot$exons$start - intron.to.zoom$end) == 
                                min(abs(transcript.to.plot$exons$start - intron.to.zoom$end)))
    
    index_second_exon <- which(abs(transcript.to.plot$exons$end - intron.to.zoom$start) == 
                                 min(abs(transcript.to.plot$exons$end - intron.to.zoom$start)))
    
  }
  
  
  exons_to_zoom <- rbind(transcript.to.plot$exons[index_first_exon,],
                         transcript.to.plot$exons[index_second_exon,])
  
  return(exons_to_zoom)
  
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



visualise_empty_plot <- function(text = NULL) {
  
  if (is.null(text)) {
    text = paste0("There are no protein-coding transcript structures available for the selected gene.")
  }
  
  ggplot() +
    theme_void() +
    geom_text(aes(x = 0, y = 0, label=text), size = 6) %>%
    return()
  
}
###############################################
## CALLS 
###############################################
# 
# create_CLIP_database_tables()
