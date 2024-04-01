
visualise_metadata_sample_clusters_gtex <- function(database_name) {
  
  # database_name <- "GTEx v8"
  
  
  db_details <- database_equivalences %>%
    filter(name == database_name)
  
  database_path <- file.path(database_folder, db_details$sqlite_file)
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  
  # Query the DB
  query <- paste0("SELECT * FROM 'metadata'")
  df_metadata <- dbGetQuery(con, query) 
  
  ## Close connection
  DBI::dbDisconnect(con)
  
  
  data_to_plot <- df_metadata %>%
    filter(cluster != "") %>%
    distinct(sample_id, .keep_all = T) %>%
    dplyr::count(cluster) %>%
    arrange(n)
  
  ggplot(data = data_to_plot, 
         aes(x=reorder(fct_rev(cluster), -n, decreasing = T), y = n, fill = cluster)) + 
    geom_bar(position = position_dodge(), stat="identity") +
    #ggforce::facet_row(~cluster, scales = "free" ) +
    #viridis::scale_fill_viridis(discrete = T, option = "F")  +
    theme_light(base_size = 16) +
    xlab("") +
    ylab("Number of samples") +
    coord_flip() +
    #ggtitle("Sample clusters and number of samples clustered within each group:")+
    #guides(fill = guide_legend(title = "Sample group: ")) +
    #geom_hline(yintercept= 70,linetype="dotted") +
    #scale_y_continuous(name ="Number of Samples", 
    #                   breaks = c(0,70,200,400,600,800))+
    
    theme(
      legend.position="none"
      #axis.text.x = element_text(color = "black",
      #                               angle = 90, 
      #                               hjust = 1,
      #                               size = "11")
    ) %>%
    return()
}

visualise_metadata_sample_clusters <- function(database_name) {
  
  # database_name <- "TCGA"
  # database_name <- "ENCODE shRNA"
  # database_name <- "LRRK2 G2019S mutation"
  
  db_details <- database_equivalences %>%
    filter(name == database_name)
  
  database_path <- file.path(database_folder, db_details$sqlite_file)
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  
  # Query the DB
  query <- paste0("SELECT * FROM 'metadata'")
  df_metadata <- dbGetQuery(con, query) 
  
  ## Close connection
  DBI::dbDisconnect(con)
  
  
  data_to_plot <- df_metadata %>%
    filter(cluster != "") %>%
    mutate(cluster = ifelse(cluster == "case", "shRNA knockdown", cluster)) %>%
    group_by(cluster) %>%
    distinct(sample_id, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(SRA_project, cluster) %>%
    arrange(n)
  data_to_plot$n %>% sum

  ggplot(data = data_to_plot, 
         aes(x=reorder(fct_rev(SRA_project), -n, decreasing = T), y = n, fill = cluster)) + 
    geom_bar(position = position_dodge(), stat="identity") +
    ggforce::facet_row(~cluster, scales = "free", 
                       labeller = label_wrap_gen(multi_line = T) ) +
    #viridis::scale_fill_viridis(discrete = T, option = "F")  +
    theme_light(base_size = 16) +
    xlab("") +
    ylab("Number of samples") +
    coord_flip() +
    #ggtitle("Sample clusters and number of samples clustered within each group:")+
    #guides(fill = guide_legend(title = "Sample group: ")) +
    #geom_hline(yintercept= 70,linetype="dotted") +
    #scale_y_continuous(name ="Number of Samples", 
    #                   breaks = c(0,70,200,400,600,800))+
    
    theme(
      legend.position="none"
      #axis.text.x = element_text(color = "black",
      #                               angle = 90, 
      #                               hjust = 1,
      #                               size = "11")
      ) %>%
    return()
}

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