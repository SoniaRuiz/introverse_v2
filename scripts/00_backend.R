library(tidyverse)
library(DBI)
library(ggplot2)
library(ggtranscript)

#setwd()

database_folder <- base::file.path("database")

database_equivalences <- data.frame(index = c(1,2,3,4,5,6),
                                    key = c("tcga",
                                            "gtex", 
                                            "encode_sh", 
                                            "encode_crispr",
                                            "pd",
                                            "ad"),
                                    name_title = c("The Cancer Genome Atlas Program (TCGA) database",
                                                   "The Genotype Tissue Expression v8 (GTEx) database",
                                                   "ENCODE shRNA database",
                                                   "ENCODE CRISPR database",
                                                   "Parkinson's Disease",
                                                   "Alzheimer's Disease"),
                                    name = c("TCGA",
                                             "GTEx v8", 
                                             "ENCODE shRNA", 
                                             "ENCODE CRISPR", 
                                             "PD/Control",
                                             "AD/Control"),
                                    sqlite_file = c("TCGA_1read_subsampleFALSE.sqlite", 
                                                    "splicing_1read.sqlite", 
                                                    "ENCODE_SR_1read_shRNA.sqlite",
                                                    "ENCODE_SR_1read_crispr.sqlite", 
                                                    "SRP058181_1read_subsampleFALSE.sqlite",
                                                    "SRP181886_1read_subsampleFALSE.sqlite"),
                                    project_types = c("cancer types",
                                                      "body sites",
                                                      "RBP knockdown experiments",
                                                      "RBP CRISPR knockdown experiments", 
                                                      "PD samples",
                                                      "AD samples"))#, "AD samples" , "PD samples" ))
                                                         

#######################################
## SET INITIAL DROPDOWN VARIABLES
#######################################


conn_list <- database_connect(database_sqlite_list = database_equivalences$sqlite_file)

set_database_project_dropdown <- function(database_key) {
  
  # message(database_key)
  # database_key <- "tcga"
  # database_key <- "gtex"
  
  database_sqlite <- database_equivalences %>% filter(key == database_key) %>% pull(sqlite_file)
  database_metadata <- get_database_metadata(database_sqlite_list = database_sqlite)[[1]] %>% as_tibble()
  
  cluster_keys <- database_metadata$SRA_project %>% unique() %>% sort()
  
  if ( database_key == "tcga" ) {
    cluster_labels <- database_metadata %>% distinct(SRA_project, tumor) %>% arrange(tumor) %>% pull(tumor)
    cluster_keys <- database_metadata %>% distinct(SRA_project, tumor) %>% arrange(tumor) %>% pull(SRA_project)
    cluster_labels <- paste0(cluster_labels, " (",cluster_keys,")")
  } else if (database_key == "encode_sh") {
    cluster_labels <- paste0("shRNA knockdown of ", cluster_keys %>% str_replace_all(pattern = "_", replacement = " "))
  } else if (database_key == "encode_crispr") {
    cluster_labels <- paste0("CRISPR editing of ", cluster_keys %>% str_replace_all(pattern = "_", replacement = " "))
  }else {
    cluster_labels <- cluster_keys %>% str_replace_all(pattern = "_", replacement = " ")  
  }
  
  
  if ( cluster_keys %>% length() > 1) {
    cluster_keys <- c("all", cluster_keys)
    cluster_labels <- c("All", cluster_labels) 
  } 
  
  names(cluster_keys) <- cluster_labels
  
  cluster_keys %>% 
    return()
  
}

set_database_sample_type_dropdown <- function(database_key, project_id) {
  
  # database_key <- "tcga"
  # project_id <- "OV"
  
  # database_key <- "encode_sh"
  # project_id <- "ADAR"
  
  # database_key <- "gtex"
  # project_id <- "BRAIN"

  message(database_key, " - ",project_id)
  
  database_sqlite <- database_equivalences %>% filter(key == database_key) %>% pull(sqlite_file)
  database_metadata <- get_database_metadata(database_sqlite_list = database_sqlite, project = project_id)[[1]] %>% as_tibble()
  
  cluster_data <- database_metadata %>%
    group_by(cluster) %>%
    distinct(sample_id, .keep_all = T) %>%
    ungroup %>%
    dplyr::count(cluster) %>%
    distinct(cluster, .keep_all = T) %>% 
    arrange(cluster) %>% 
    drop_na() %>%
    filter(cluster != "") %>%
    as.data.frame() %>%
    mutate(cluster_label = paste0(cluster %>% str_to_sentence(), " (",n," samples)"))
  
  if ( cluster_data %>% nrow() > 1 ) {
    cluster_data <- rbind(data.frame(cluster = "all", n = "", cluster_label = "All"), cluster_data)
  } 
  
  if ( any(str_detect(string = str_to_lower(cluster_data$cluster_label), pattern = "case")) ) {
    index_case <- which(str_detect(string = str_to_lower(cluster_data$cluster_label),
                                   pattern = "case"))
    
    if (project_id == "encode_sh") {
      cluster_label = "shRNA knockdown"
    } else {
      cluster_label = "CRISPR editing"
    }
    cluster_data$cluster_label[index_case] <- str_replace(string = cluster_data$cluster_label[index_case] %>% str_to_lower(),
                                                          pattern = "case",
                                                          replacement = cluster_label)
  }
  
  cluster_keys <- cluster_data$cluster
  names(cluster_keys) <- cluster_data$cluster_label

  
  cluster_keys %>% return()
}

## Gene List
gene_list <- get_database_genes(conn_list)

