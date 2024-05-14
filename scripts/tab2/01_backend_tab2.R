


###################################################
## GENE QUERY - SPLICING VISUALISATION
###################################################


get_gene_splicing_data_to_visualise <- function(gene.id, 
                                                transcript.id,
                                                database.sqlite, 
                                                junction.type, 
                                                project.id = NULL,
                                                table.name = NULL,
                                                compare = NULL) {
  
  
  
  message(compare, " - ", gene.id, " - ", transcript.id, " - ", database.sqlite, " - ", project.id, " - ", table.name, " - ", junction.type)
  
  
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
  
  # compare = ""
  # gene.id = "PTEN"
  # transcript.id = "ENST00000371953"
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id = "CESC"
  # table.name = "Primary tumor"
  # junction.type = "all"
  
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
  DBI::dbListTables(con)
  #DBI::dbGetQuery(con,  paste0("SELECT * FROM 'intron'")) %>% as_tibble()
  
  
  #####################################################
  ## BUILD SPLICING QUERY 
  #####################################################
  
  ## Merge gene and transcript info
  query <- paste0("SELECT transcript.id, transcript.transcript_id AS transcript_ENS, transcript.MANE, gene.gene_name 
                  FROM 'transcript' INNER JOIN 'gene' ON transcript.gene_id=gene.id
                  WHERE gene.gene_name = '", gene.id, "'")
  db_gene_w_transcript <- DBI::dbGetQuery(con, query) %>% as_tibble()
  
  
  
  ## Get mis-splicing from selected sample group
  query <- paste0("SELECT missplicing.ref_junID, missplicing.novel_junID, missplicing.transcript_id, missplicing.novel_sum_counts, missplicing.ref_sum_counts,
                  novel.novel_type, novel.seqnames AS novel_seqnames, novel.start AS novel_start, novel.end AS novel_end, novel.strand AS novel_strand,
                  intron.seqnames AS intron_seqnames, intron.start AS intron_start, intron.end AS intron_end, intron.strand AS intron_strand,
                  intron.transcript_id_list
                  FROM '", table.name, "_", project.id, "_misspliced' AS 'missplicing'
                  INNER JOIN 'novel' ON novel.novel_junID = missplicing.novel_junID
                  INNER JOIN 'intron' ON intron.ref_junID = missplicing.ref_junID
                  WHERE missplicing.transcript_id IN ('",paste(db_gene_w_transcript$id, collapse="','"),"' )")
  db_sample_cluster_missplicing <- DBI::dbGetQuery(con, query) %>% as_tibble()
  
  
  
  ## Get accurate splicing from selected sample group
  query <- paste0("SELECT nevermisspliced.ref_junID, nevermisspliced.transcript_id, nevermisspliced.ref_sum_counts,
                  intron.seqnames AS intron_seqnames, intron.start AS intron_start, intron.end AS intron_end, intron.strand AS intron_strand,
                  intron.transcript_id_list
                  FROM '", table.name, "_", project.id, "_nevermisspliced' AS nevermisspliced
                  INNER JOIN 'intron' ON intron.ref_junID = nevermisspliced.ref_junID
                  WHERE nevermisspliced.transcript_id IN ('",paste(db_gene_w_transcript$id, collapse="','"),"')")
  db_sample_cluster_splicing <- DBI::dbGetQuery(con, query) %>% as_tibble()
  
  
  
  ## Merge mis-splicing and accurate splicing
  db_sample_cluster_all_splicing = plyr::rbind.fill(db_sample_cluster_missplicing, db_sample_cluster_splicing) %>% as_tibble()  %>%
    dplyr::select(-transcript_id) %>%
    mutate(transcript_id = str_split(transcript_id_list, pattern = ",")) %>%
    unnest(transcript_id) %>%
    mutate(transcript_id = transcript_id %>% as.integer()) %>%
    inner_join(y = db_gene_w_transcript,
               by = c("transcript_id" = "id")) %>%
    filter(transcript_ENS == transcript.id )
  
  if ( !is.null(compare) && compare != "") {
    db_sample_cluster_all_splicing <- db_sample_cluster_all_splicing %>%
      filter(novel_junID %in% compare)
  }
  
  if ( db_sample_cluster_all_splicing %>% nrow() == 0 ) {
    return(list(error = paste0("Sorry, there is no splicing activity associated with '", gene.id, "' in the sample group and experiment selected.")))
  } else {
    #####################################################
    ## VISUALISE
    #####################################################
    
    ## Get coordinates from the transcript selected
    hg38_transcripts_to_plot <- get_transcript_to_plot(transcript_id = transcript.id,
                                                       multiple = T)  
    
    if ( is.null(hg38_transcripts_to_plot) | 
         (!is.null(hg38_transcripts_to_plot) && (nrow(hg38_transcripts_to_plot$cds) == 0 | nrow(hg38_transcripts_to_plot$utr) == 0)) ) {
      
      # print(hg38_transcripts_to_plot$exons$transcript_biotype %>% unique)
      return(list(error = paste0("Sorry, it has not been possible to visualise the selected transcript structure.\nIt might not have CDS or UTR regions defined.\nTranscript biotype (Ensembl v111): '", 
                                 str_replace(string = hg38_transcripts_to_plot$exons$transcript_biotype %>% unique, pattern = "_", replacement = " "), "'.")))
      
      
    } else {
      
      ## PLOT
      
      if ( hg38_transcripts_to_plot$exons %>% nrow() > 0 ) {
        
        
        annotated_introns_to_plot <- db_sample_cluster_all_splicing %>%
          dplyr::select(seqnames = intron_seqnames,
                        start = intron_start,
                        end = intron_end,
                        strand = intron_strand,
                        transcript_id = transcript_ENS,
                        
                        sum_counts = ref_sum_counts) %>%
          distinct(start, end, transcript_id, .keep_all = T) %>%
          mutate(junction_color = "blue",
                 junction_type = "annotated_intron") %>%
          as_tibble() %>%
          filter(transcript_id %in% (hg38_transcripts_to_plot$exons$transcript_id %>% unique)) %>%
          filter(start %in% (ggtranscript::to_intron(hg38_transcripts_to_plot$exons, "transcript_id")$start +1 ), 
                 end %in% (ggtranscript::to_intron(hg38_transcripts_to_plot$exons, "transcript_id")$end -1 ))
        
        
        novel_junctions_to_plot <- db_sample_cluster_all_splicing %>% 
          dplyr::select(seqnames = novel_seqnames, start = novel_start, end = novel_end, strand = novel_strand, 
                        transcript_id = transcript_ENS, 
                        junction_type = novel_type, 
                        sum_counts = novel_sum_counts) %>%
          #filter(novel_sum_counts > 1) %>%
          distinct(start, end, transcript_id, .keep_all = T) %>% 
          mutate(length = end - start,
                 junction_color = ifelse(junction_type == "novel_donor", "#288a5b", "#64037d")) %>%
          as_tibble() %>%
          filter(length < sum(hg38_transcripts_to_plot$transcript$width))%>%
          dplyr::select(-length) %>%
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
        
        
        
        if ( junctions_to_plot %>% nrow > 0 ) {
          
          list(exons = hg38_transcripts_to_plot$exons,
               cds = hg38_transcripts_to_plot$cds,
               utr = hg38_transcripts_to_plot$utr,
               junctions_to_plot = junctions_to_plot) %>%
            return()
          
          
        } else {
          list(error = paste0("Sorry, there are no splicing events meeting the selected criteria.")) %>% return()
        }
        
        
        
      } else {
        list(error = paste0("Sorry, there are no splicing events meeting the selected criteria.")) %>% return()
      }
    }
  }
  
  
  
  #})
  #query_results %>% return()
}




get_required_data_before_gene_comparison <- function(database.sqlite,
                                                     project.id,
                                                     sample.type, 
                                                     compare.table,
                                                     junction.type,
                                                     gene.name) {
  
  
  # database.sqlite = "ENCODE_SR_1read.sqlite"
  # project.id ="ADAR"
  # sample.type = "all" 
  # compare.table = "all"
  # gene.name = "APOE"
  
  # gene.name = "MAPT"
  
  
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id ="BLCA"
  # sample.type = "Primary tumor"
  # compare.table = "all"
  # gene.name = "A1BG-AS1"
  
  
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id ="ACC"
  # sample.type = "Primary tumor"
  # compare.table = "all"
  # gene.name = "A1BG-AS1"
  
  # database.sqlite = "TCGA_1read_subsampleFALSE.sqlite"
  # project.id ="BLCA"
  # sample.type = "Primary tumor"
  # compare.table = "all"
  # gene.name = "PTEN"
  
  message(  database.sqlite, " ", project.id, " ", sample.type, " ", compare.table, " ", gene.name, " ")
  
  
  unique_junctions_to_process = NULL
  unique_splicing <- list()
  
  if ( compare.table == "compare" ){
    
    
    database_table_names <- (get_database_metadata(database_sqlite_list  = database.sqlite, 
                                                   project = project.id))[[1]] %>% 
      dplyr::distinct(cluster) %>% 
      filter(cluster != "") %>% 
      pull()
    
    
    for (table in database_table_names) {
      
      # table <- database_table_names[1]
      query <- paste0("SELECT misspliced.ref_junID, misspliced.novel_junID
                      FROM '", table, "_", project.id, "_misspliced' AS misspliced
                      INNER JOIN transcript ON transcript.id=misspliced.transcript_id
                      INNER JOIN gene ON gene.id = transcript.gene_id
                      WHERE gene_name ='", gene.name, "'")
      db_novel <- DBI::dbGetQuery(conn_list[database.sqlite][[1]], query) %>% as_tibble()
      
      
      # table <- database_table_names[1]
      query <- paste0("SELECT never.ref_junID
                      FROM '", table, "_", project.id, "_nevermisspliced' AS never
                      INNER JOIN transcript ON transcript.id=never.transcript_id
                      INNER JOIN gene ON gene.id = transcript.gene_id
                      WHERE gene_name ='", gene.name, "'")
      db_never <- DBI::dbGetQuery(conn_list[database.sqlite][[1]], query) %>% as_tibble()
      
      
      if (junction.type == "all") {
        unique_splicing[[table]] = c(db_novel$novel_junID, db_novel$ref_junID, db_never$ref_junID)
      } else if (junction.type == "intron") {
        unique_splicing[[table]] = c(db_novel$ref_junID, db_never$ref_junID)
      } else {
        unique_splicing[[table]] = db_novel$novel_junID
      }
      
      
    }
    
    ## Get index of table to perform unique junx
    ind <- which(names(unique_splicing) == sample.type)
    
    if ( length(names(unique_splicing[-c(ind)])) > 0 ) {
      unique_junctions_to_process <- setdiff(unique_splicing[ind][[1]],
                                             do.call(c, unique_splicing[-ind]) %>% unname())
    } else {
      unique_junctions_to_process <- setdiff(unique_splicing[[ind]],
                                             unique_splicing[[-ind]])
    }

  } 
  
  return(list(unique_junctions = unique_junctions_to_process,
              table_name = sample.type,
              compare = ifelse(compare.table == "compare", T, F)))
}


visualise_gene_data <- function(gene_splicing_data) {
  
  
  if (any(names(gene_splicing_data) == "error")) {
    
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
      # ggtranscript::geom_junction(
      #   data = annotated_introns_to_plot,
      #   color = "blue",
      #   ncp = 100,
      #   junction.y.max = 0.5
      # ) +
      ggtranscript::geom_junction(
        data = gene_splicing_data$junctions_to_plot,
        color = gene_splicing_data$junctions_to_plot$junction_color,
        ncp = 100,
        junction.y.max = 0.5
      ) +
      ggtranscript::geom_junction_label_repel(
        data = gene_splicing_data$junctions_to_plot,
        mapping = aes(label = sum_counts),
        
        max.overlaps = 150,
        junction.y.max = 0.5
      ) +
      
      theme_light(base_size = 16) +
      ggtitle(ifelse(table.name == "case", "shRNA knockdown", table.name))+
      
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
            legend.position = "top",
            legend.text = element_text(size = "11"),
            legend.title = element_text(size = "11")) +
      xlab(paste0("Genomic position (chr", junctions_to_plot$exons$seqnames %>% unique ,")")) + 
      ylab("") %>% 
      return()
  }
}