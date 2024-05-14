



# 1. Visualise number of samples per cluster
visualise_metadata_sample_clusters_gtex <- function(data_to_plot) {
  
  # database_name <- "GTEx v8"
  
  ggplot(data = data_to_plot, 
         aes(x = reorder(fct_rev(cluster), -n, decreasing = T), 
             y = n, 
             fill = SRA_project)) + 
    geom_bar(position = position_dodge(), stat="identity") +
    theme_light(base_size = 16) +
    xlab("") +
    ylab("Number of samples") +
    coord_flip() +
    
    theme(
      legend.position="none"
    ) %>%
    return()
}

visualise_metadata_sample_clusters_encode <- function(data_to_plot){
  
  
  ggplot(data = data_to_plot, 
         aes(x=reorder(SRA_project, -n, decreasing = T), y = n, fill = cluster)) + 
    geom_bar(position = position_dodge(), stat="identity") +
    ggforce::facet_col(~cluster,  
                       labeller = labeller(cluster = c("control" = "Control",
                                                       "case" = "shRNA knockdown"))) +
    theme_light(base_size = 16) +
    xlab("") +
    ylab("Number of samples") +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1),
                                                                n=2))))+
    theme(
      legend.position="none",
      axis.text.x = element_text(color = "black",
                                     angle = 90, 
                                     hjust = 1, vjust = 0.5,
                                   size = "11")
    ) %>%
    return()
}

visualise_metadata_sample_clusters_SRP100948 <- function(data_to_plot){
  
  
  ggplot(data = data_to_plot, 
         aes(x=reorder(cluster, -n, decreasing = T), y = n, fill = cluster)) + 
    geom_bar(position = position_dodge(), stat="identity") +
    ggforce::facet_row(~SRA_project) +
    theme_light(base_size = 16) +
    xlab("") + 
    ylab("Number of samples") +
    theme(legend.position="none") +
    scale_y_continuous(n.breaks = 6) %>%
    return()
}

visualise_metadata_sample_clusters_SRP151040 <- function(data_to_plot){
  
  
  ggplot(data = data_to_plot, 
         aes(x=reorder(cluster, -n, decreasing = T), y = n, fill = cluster)) + 
    geom_bar(position = position_dodge(), stat="identity") +
    ggforce::facet_row(~SRA_project) +
    theme_light(base_size = 16) +
    xlab("") + 
    ylab("Number of samples") +
    theme(legend.position="none") +
    coord_flip() %>%
    return()
}

visualise_metadata_sample_clusters <- function(database_name) {
  
  # database_name <- "TCGA"
  # database_name <- "GTEx v8"
  # database_name <- "GTEx v8 Age"
  # database_name <- "ENCODE shRNA"
  # database_name <- "LRRK2 G2019S mutation"
  # database_name <- "Alzheimer's Disease/Control"
  # database_name <- "Parkinson's Disease/Control"
  # database_name <- "LRRK2 G2019S mutation" 
  
  # ####################################
  # ## GET THE DATABASE METADATA
  # ####################################

  db_details <- database_equivalences %>%
    filter(name == database_name)

  # Connect to the database
  database_path <- file.path(database_folder, db_details$sqlite_file)
  con <- dbConnect(RSQLite::SQLite(), database_path)

  # Query the DB
  query <- paste0("SELECT * FROM 'metadata'")
  df_metadata <- dbGetQuery(con, query)

  ## Close connection
  DBI::dbDisconnect(con)
  
  
  ####################################
  ## TIDY DATA BEFORE PLOT
  ####################################
  
  data_to_plot <- df_metadata %>%
    filter(cluster != "") %>%
    group_by(SRA_project, cluster) %>%
    distinct(sample_id, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(cluster, SRA_project, sort = TRUE) %>%
    mutate(SRA_project = stringr::str_replace_all(string = SRA_project, pattern = "_", replacement = " "),
           SRA_project = SRA_project %>% as.factor(),
           cluster = cluster %>% as.factor()) 

  
  ####################################
  ## PLOT
  ####################################
  
  if ( database_name == "GTEx v8") {
    visualise_metadata_sample_clusters_gtex(data_to_plot)
  } else if ( database_name == "ENCODE shRNA") {
    visualise_metadata_sample_clusters_encode(data_to_plot)
  } else if ( database_name == "Alzheimer's Disease/Control" ||
              database_name == "Parkinson's Disease/Control" ) {
    visualise_metadata_sample_clusters_SRP100948(data_to_plot)
  } else if ( database_name == "LRRK2 G2019S mutation"  ) {
    visualise_metadata_sample_clusters_SRP151040(data_to_plot)
  } else {
    
    data_to_plot <- data_to_plot %>%
      mutate(SRA_project = tidytext::reorder_within(SRA_project, n, cluster)) 
    
    data_to_plot %>%
      ggplot( aes(y = SRA_project,  x = n, fill = cluster) ) +
      geom_bar(stat = 'identity', show.legend = FALSE) +
      ggforce::facet_row(~cluster, scales = "free", labeller = label_wrap_gen()) +
      tidytext::scale_y_reordered() +
      theme_light(base_size = 16) +
      ylab("") +
      xlab("Number of samples") +
      #coord_flip() +
      #scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1), n=2))))+
      theme(legend.position="none") %>%
      return()
  }
  
}


# 2. Visualise age per cluster

table_cancer_types <- function() {
  
  data.frame("Study Abbreviation" = c("LAML","ACC ","BLCA","LGG ","BRCA","CESC","CHOL","LCML","COAD","ESCA","FPPP","GBM ","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","DLBC","MESO","MISC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THYM","THCA","UCS ","UCEC","UVM"),
             "Study Name" = c("Acute Myeloid Leukemia","Adrenocortical carcinoma","Bladder Urothelial Carcinoma","Brain Lower Grade Glioma","Breast invasive carcinoma","Cervical squamous cell carcinoma and endocervical adenocarcinoma",
                              "Cholangiocarcinoma","Chronic Myelogenous Leukemia","Colon adenocarcinoma","Esophageal carcinoma","FFPE Pilot Phase II","Glioblastoma multiforme","Head and Neck squamous cell carcinoma",
                              "Kidney Chromophobe","Kidney renal clear cell carcinoma","Kidney renal papillary cell carcinoma","Liver hepatocellular carcinoma","Lung adenocarcinoma","Lung squamous cell carcinoma",
                              "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma","Mesothelioma","Miscellaneous","Ovarian serous cystadenocarcinoma","Pancreatic adenocarcinoma","Pheochromocytoma and Paraganglioma",
                              "Prostate adenocarcinoma","Rectum adenocarcinoma","Sarcoma","Skin Cutaneous Melanoma","Stomach adenocarcinoma","Testicular Germ Cell Tumors","Thymoma","Thyroid carcinoma",
                              "Uterine Carcinosarcoma","Uterine Corpus Endometrial Carcinoma","Uveal Melanoma")) %>%
    return()
}

