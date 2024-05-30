


###################################################
## GENE QUERY - SPLICING VISUALISATION
###################################################


get_unique_jxn_from_comparison <- function(gene.id,
                                           transcript.id,
                                           database.sqlite,
                                           project.id,
                                           table.name, 
                                           junction.type) {
  
  
  # database.sqlite = "ENCODE_SR_1read.sqlite"
  # project.id ="ADAR"
  # table.name = "all" 
  # gene.id = "APOE"
  
  # gene.id = "MAPT"
  
  
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id ="BLCA"
  # table.name = "Primary tumor"
  # gene.id = "A1BG-AS1"
  
  
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id ="ACC"
  # table.name = "Primary tumor"
  # gene.id = "A1BG-AS1"
  
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id ="BLCA"
  # table.name = "Primary tumor"
  # gene.id = "PTEN"
  
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id ="SARC"
  # table.name = "Primary tumor"
  # table.name = "Solid tissue normal"
  # gene.id = "A2M"
  # junction.type = "novel_acceptor"
  # junction.type = "intron"
  # transcript.id = "ENST00000318602"
  
  
  # gene.id = "GBA1"
  # transcript.id = "ENST00000368373"
  # database.sqlite = "splicing_1read.sqlite"
  # project.id = "BRAIN"
  # table.name = "Brain - Cortex"
  # junction.type = "all"
  message(  gene.id, " | ", transcript.id, " | ", database.sqlite, " | ", project.id, " | ", table.name, " | ", junction.type)
  
  
  unique_junctions_to_process = NULL
  unique_splicing <- list()
  

  database_table_names <- (get_database_metadata(database_sqlite_list  = database.sqlite, project = project.id))[[1]] %>% 
    dplyr::distinct(cluster) %>% 
    filter(cluster != "") %>% 
    pull()
  
  
  for (table in database_table_names) {
    
    # table <- database_table_names[1]
    # table <- database_table_names[2]
    
    message(table, "...")
    query <- paste0("SELECT misspliced.ref_junID, misspliced.novel_junID
                    FROM '", table, "_", project.id, "_misspliced' AS misspliced
                    INNER JOIN transcript ON transcript.id=misspliced.transcript_id
                    INNER JOIN gene ON gene.id = transcript.gene_id
                    WHERE gene_name ='", gene.id, "' AND transcript.transcript_id='", transcript.id, "'")
    db_novel <- DBI::dbGetQuery(conn_list[database.sqlite][[1]], query) %>% as_tibble()
    
    
    # table <- database_table_names[1]
    query <- paste0("SELECT never.ref_junID
                    FROM '", table, "_", project.id, "_nevermisspliced' AS never
                    INNER JOIN transcript ON transcript.id=never.transcript_id
                    INNER JOIN gene ON gene.id = transcript.gene_id
                    WHERE gene_name ='", gene.id, "' AND transcript.transcript_id='", transcript.id, "'")
    db_never <- DBI::dbGetQuery(conn_list[database.sqlite][[1]], query) %>% as_tibble()
    
    
    if (junction.type == "all") {
      unique_splicing[[table]] = c(db_novel$novel_junID, db_novel$ref_junID, db_never$ref_junID) %>% unique()
    } else if (junction.type == "intron") {
      unique_splicing[[table]] = c(db_novel$ref_junID, db_never$ref_junID) %>% unique()
    } else {
      unique_splicing[[table]] = db_novel$novel_junID %>% unique()
    }
    
  }
  
  
  ## Get index of table to perform unique junx
  ind <- which(names(unique_splicing) == table.name)
  
  if ( length(names(unique_splicing[-c(ind)])) > 0 ) {
    unique_junctions_to_process <- setdiff(unique_splicing[ind][[1]] %>% unique,
                                           (do.call(c, unique_splicing[-ind]) %>% unname() %>% unique))
  } else {
    unique_junctions_to_process <- setdiff(unique_splicing[[ind]],
                                           unique_splicing[[-ind]])
  }
    

  
  return(list(unique_junctions = unique_junctions_to_process,
              table_name = table.name))
}


# compare = F
# gene.id = "A2M"
# transcript.id = "ENST00000318602"
# database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
# project.id = "LAML"
# table.name = "Primary blood derived cancer - peripheral blood"
# junction.type = "all"
# 
# 
# gene_data_to_visualise <- get_gene_splicing_data_to_visualise(gene.id,
#                                                               transcript.id,
#                                                               database.sqlite,
#                                                               junction.type,
#                                                               project.id,
#                                                               table.name,
#                                                               compare)


get_gene_splicing_data_to_visualise <- function(gene.id, 
                                                transcript.id,
                                                database.sqlite, 
                                                project.id,
                                                table.name, 
                                                junction.type,
                                                compare) {
  
  
  
  message(compare, " - ", gene.id, " - ", transcript.id, " - ", database.sqlite, " - ", project.id, " - ", table.name, " - ", junction.type)
  
  # gene.id = "PSEN1"
  # transcript.id = "ENST00000394164"
  # database.sqlite="ENCODE_SR_1read.sqlite"
  # project.id = "TARDBP"
  # table.name = "case"
  # junction.type="all"
  # compare = F
  
  # gene.id="PTER"
  # transcript.id="ENST00000535784"
  # database.sqlite="ENCODE_SR_1read.sqlite"
  # project.id="TARDBP"
  # table.name="case"
  # junction.type="all"
  # compare = F
  
  # gene.id="HDGFL2"
  # transcript.id="ENST00000616600"
  # database.sqlite="ENCODE_SR_1read.sqlite"
  # project.id="TARDBP"
  # table.name="case"
  # junction.type="all"
  # compare = T
  
  # compare = ""
  # gene.id = "MAPT"
  # transcript.id = "ENST00000371953"
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id = "LGG"
  # table.name = "Primary tumor"
  # junction.type = "all"
  
  # compare = ""
  # gene.id = "PTEN"
  # transcript.id = "ENST00000371953"
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id = "ACC"
  # table.name = "Primary tumor"
  # junction.type = "all"
  
  # compare = F
  # gene.id = "PTEN"
  # transcript.id = "ENST00000371953"
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id = "CESC"
  # table.name = "Primary tumor"
  # junction.type = "intron"
  
  # gene.id = "APOE"
  # transcript.id = "ENST00000252486"
  
  
  # gene.id = "PSEN1"
  # transcript.id = "ENST00000324501"
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id = "ACC"
  # table.name = "Primary tumor"
  
  # gene.id = "PTEN" 
  # database.sqlite = "ENCODE_SR_1read.sqlite"
  # project.id = "TARDBP"
  # table.name = "case"
  # compare = c(97656, 109693)
  
  
  # gene.id = "HDGFL2"
  # transcript.id = "ENST00000616600"
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id = "CESC"
  # table.name = "Solid tissue normal"
  # compare = 1370779
  
  
  # gene.id = "A1BG-AS1"
  # transcript.id = "ENST00000593374"
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id = "BLCA"
  # table.name = "Primary tumor"
  # compare = ""
  # print(databases_name)
  
  # if ( any(databases_name == "all") ){
  #   databases_sqlite_list <- database_equivalences$sqlite_file
  #   project.id <- NULL
  #   table.name <- NULL
  # } else {
  #   databases_sqlite_list <- database_equivalences %>% filter(key == databases_name) %>% pull(sqlite_file)
  # }
  
  
  
  con <- conn_list[database.sqlite][[1]]
  #DBI::dbListTables(con) %>% print
  
  if (compare) {
    
    unique_jxn_from_comparison <- get_unique_jxn_from_comparison(gene.id = gene.id,
                                                                 transcript.id = transcript.id,
                                                                 database.sqlite = database.sqlite,
                                                                 project.id = project.id,
                                                                 table.name = table.name, 
                                                                 junction.type = junction.type)
    
    if ( length(unique_jxn_from_comparison$unique_junctions) == 0 ) {
      
      return(list(error = paste0("There are no ", 
                                 ifelse(junction.type == "all", "",
                                        str_replace_all(string =  paste0(junction.type, " "),
                                                        pattern = "_", 
                                                        replacement = " ")),
                                 "splicing events originating from the '", gene.id, 
                                 "' gene\nthat are only found across '", table.name , 
                                 "' samples\nwithin the '", project.id,"' experiment.")))
    }
  }
  
  
  
  
  #####################################################
  ## BUILD SPLICING QUERY 
  #####################################################
  
  ## Merge gene and transcript info
  query <- paste0("SELECT transcript.id, transcript.transcript_id AS transcript_ENS, transcript.MANE, gene.gene_name 
                  FROM 'transcript' INNER JOIN 'gene' ON transcript.gene_id=gene.id
                  WHERE gene.gene_name = '", gene.id,"'")
  db_gene_w_transcript <- DBI::dbGetQuery(con, query) %>% as_tibble()
  
  
  
  ## Get mis-splicing from selected sample group
  query <- paste0("SELECT missplicing.ref_junID, missplicing.novel_junID, missplicing.transcript_id, missplicing.novel_sum_counts, missplicing.ref_sum_counts,
                  missplicing.novel_n_individuals, missplicing.ref_n_individuals,
                  novel.novel_type, novel.seqnames AS novel_seqnames, novel.start AS novel_start, novel.end AS novel_end, novel.strand AS novel_strand,
                  intron.seqnames AS intron_seqnames, intron.start AS intron_start, intron.end AS intron_end, intron.strand AS intron_strand,
                  intron.transcript_id_list
                  FROM '", table.name, "_", project.id, "_misspliced' AS 'missplicing'
                  INNER JOIN 'novel' ON novel.novel_junID = missplicing.novel_junID
                  INNER JOIN 'intron' ON intron.ref_junID = missplicing.ref_junID
                  WHERE missplicing.transcript_id IN ('",paste(db_gene_w_transcript$id, collapse="','"),"' )")
  db_sample_cluster_missplicing <- DBI::dbGetQuery(con, query) %>% as_tibble()
  
  
  
  ## Get accurate splicing from selected sample group
  query <- paste0("SELECT nevermisspliced.ref_junID, nevermisspliced.transcript_id, nevermisspliced.ref_sum_counts, nevermisspliced.ref_n_individuals,
                  intron.seqnames AS intron_seqnames, intron.start AS intron_start, intron.end AS intron_end, intron.strand AS intron_strand,
                  intron.transcript_id_list
                  FROM '", table.name, "_", project.id, "_nevermisspliced' AS nevermisspliced
                  INNER JOIN 'intron' ON intron.ref_junID = nevermisspliced.ref_junID
                  WHERE nevermisspliced.transcript_id IN ('",paste(db_gene_w_transcript$id, collapse="','"),"')")
  db_sample_cluster_nevermisspliced <- DBI::dbGetQuery(con, query) %>% as_tibble()
  

  
  ## Merge mis-splicing and accurate splicing
  db_sample_cluster_all_splicing_w_transcript <- plyr::rbind.fill(db_sample_cluster_missplicing, db_sample_cluster_nevermisspliced) %>% 
    as_tibble()  %>%
    dplyr::select(-transcript_id) %>%
    mutate(transcript_id = str_split(transcript_id_list, pattern = ",")) %>%
    unnest(transcript_id) %>%
    mutate(transcript_id = transcript_id %>% as.integer()) %>%
    inner_join(y = db_gene_w_transcript, by = c("transcript_id" = "id")) %>%
    filter(transcript_ENS == transcript.id )
  
  
  
 
  
  if ( db_sample_cluster_all_splicing_w_transcript %>% nrow() == 0 ) {
    
    return(list(error = paste0("Sorry, there is no splicing activity associated with '", gene.id, "' in the sample group and experiment selected.")))
    
  } else {
    
    #####################################################
    ## VISUALISE
    #####################################################
    
    ## Get coordinates from the transcript selected
    
    hg38_transcripts_to_plot <- get_transcript_to_plot(transcript.id = transcript.id,
                                                       jxn.type = junction.type,
                                                       multiple = T)  
    
    ## PLOT
    
    if ( !is.null(hg38_transcripts_to_plot) && hg38_transcripts_to_plot$exons %>% nrow() > 0 ) {
      
      
      annotated_introns_to_plot <- db_sample_cluster_all_splicing_w_transcript %>%
        dplyr::select(ID = ref_junID,
                      seqnames = intron_seqnames,
                      start = intron_start,
                      end = intron_end,
                      strand = intron_strand,
                      transcript_id = transcript_ENS,
                      sum_counts = ref_sum_counts,
                      "Seen in (n samples)" = ref_n_individuals) %>%
        arrange(start) %>%
        distinct(start, end, transcript_id, .keep_all = T) %>%
        mutate(junction_color = "blue",
               junction_type = "annotated_intron") %>%
        as_tibble() %>%
        filter(transcript_id %in% (hg38_transcripts_to_plot$exons$transcript_id %>% unique)) %>%
        filter(start %in% (ggtranscript::to_intron(hg38_transcripts_to_plot$exons, "transcript_id")$start + 1 ), 
               end %in% (ggtranscript::to_intron(hg38_transcripts_to_plot$exons, "transcript_id")$end - 1 )) 
      

      
      novel_junctions_to_plot <- db_sample_cluster_all_splicing_w_transcript %>% 
        dplyr::select(seqnames = novel_seqnames, start = novel_start, end = novel_end, strand = novel_strand, 
                      ID = novel_junID,
                      transcript_id = transcript_ENS, 
                      junction_type = novel_type, 
                      sum_counts = novel_sum_counts,
                      "Seen in (n samples)" = novel_n_individuals) %>%
        #filter(novel_sum_counts > 1) %>%
        distinct(start, end, transcript_id, .keep_all = T) %>% 
        mutate(length = end - start,
               junction_color = ifelse(junction_type == "novel_donor", "#288a5b", "#64037d")) %>%
        as_tibble() %>%
        filter(length < sum(hg38_transcripts_to_plot$transcript$width))%>%
        dplyr::select(-length) %>%
        arrange(start) %>%
        filter(transcript_id %in% (hg38_transcripts_to_plot$exons$transcript_id %>% unique)) %>%
        filter(start %in% (ggtranscript::to_intron(hg38_transcripts_to_plot$exons, "transcript_id")$start +1 ) |
                 end %in% (ggtranscript::to_intron(hg38_transcripts_to_plot$exons, "transcript_id")$end -1 ))
      
      
      if (junction.type == "all") {
        junctions_to_plot <- plyr::rbind.fill(annotated_introns_to_plot, novel_junctions_to_plot)
      } else if (junction.type == "intron") {
        junctions_to_plot <- annotated_introns_to_plot
      } else {
        junctions_to_plot <- novel_junctions_to_plot %>% filter(junction_type == junction.type)
      }
      
      # ggtranscript::to_intron(hg38_transcripts_to_plot$exons, "transcript_id")%>% arrange(start + 1)
      # junctions_to_plot%>% arrange(start)
      
      if ( compare ) {
        
        junctions_to_plot <- junctions_to_plot %>%
          filter(ID %in% unique_jxn_from_comparison$unique_junctions)
        
      }

      
      if ( nrow(hg38_transcripts_to_plot$cds) == 0 | nrow(hg38_transcripts_to_plot$utr) == 0 ) {
        
        error_message <- paste0("Sorry, it has not been possible to visualise the selected transcript structure.\nIt might not have CDS or UTR regions defined.\nTranscript biotype (Ensembl v111): '", 
                                str_replace(string = hg38_transcripts_to_plot$exons$transcript_biotype %>% unique, pattern = "_", replacement = " "), "'.")
        
        if ( junctions_to_plot %>% nrow > 0 ) {
          return(list(error = error_message,
                      junctions_to_plot = junctions_to_plot,
                      plot_title = table.name))
        } else {
          return(list(error = error_message,
                      junctions_to_plot = NULL,
                      plot_title = NULL))
        }
        
        
      } else {
        
        if ( junctions_to_plot %>% nrow > 0 ) {
          
          list(exons = hg38_transcripts_to_plot$exons,
               cds = hg38_transcripts_to_plot$cds,
               utr = hg38_transcripts_to_plot$utr,
               junctions_to_plot = junctions_to_plot,
               plot_title = table.name) %>%
            return()
          
          
        } else {
          
          list(error = paste0("There are no ", 
                              ifelse(junction.type == "all", "",
                                     str_replace_all(string =  paste0(junction.type, " "),
                                                     pattern = "_", 
                                                     replacement = " ")), 
                              "splicing events originating from the '", 
                              gene.id, "' gene\nuniquely found across '", table.name , 
                              "' samples\nwithin the '", project.id,"' experiment.")) %>% return()
        }
        
      }
      
      
      
    } else {
      list(error = paste0("Sorry, there are no splicing events meeting the selected criteria.")) %>% return()
    }
  }
  
  
  
  
  #})
  #query_results %>% return()
}


visualise_gene_data <- function(gene_splicing_data) {
  
  
  if ( any(names(gene_splicing_data) == "error") ) {
    
    visualise_empty_plot(gene_splicing_data$error) %>%
      return()
    
  } else {
    
    gene_splicing_data$exons %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = transcript_id 
      )) +
      ggtranscript::geom_range(
        data = gene_splicing_data$cds,
        fill = "#999999"
      ) +
      ggtranscript::geom_range(
        data = gene_splicing_data$utr,
        fill = "#333333",
        height = 0.25
      ) +
      ggtranscript::geom_intron(
        data = ggtranscript::to_intron(gene_splicing_data$exons, "transcript_id"),
        aes(strand = strand %>% unique)
      ) + 
      ggtranscript::geom_junction(
        data = gene_splicing_data$junctions_to_plot,
        color = gene_splicing_data$junctions_to_plot$junction_color,
        ncp = 100, 
        junction.y.max = 0.5
      ) +
      ggtranscript::geom_junction_label_repel(
        data = gene_splicing_data$junctions_to_plot,
        color = gene_splicing_data$junctions_to_plot$junction_color,
        mapping = aes(label = sum_counts), 
        ncp = 100, 
        junction.y.max = 0.5
        
      ) + 
      
      theme_light(base_size = 16) +
      ggtitle(ifelse(gene_splicing_data$plot_title == "case", 
                     "shRNA knockdown", gene_splicing_data$plot_title))+
      
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
            legend.position = "top",
            legend.text = element_text(size = "11"),
            legend.title = element_text(size = "11")) +
      xlab(paste0("Genomic position (", gene_splicing_data$exons$seqnames %>% unique ,")")) + 
      ylab(gene_splicing_data$exons$gene_name %>% unique) %>% 
      return()
  }
}

