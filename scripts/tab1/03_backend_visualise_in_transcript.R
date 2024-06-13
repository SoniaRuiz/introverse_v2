
visualise_transcript <- function(junID,
                                 transcript_ENS,
                                 database_name,
                                 junType) {


  # junType <- 'Novel Donor'
  # junID = "chr1:155235312-155235462:-"
  # transcript_ENS = "ENST00000368373"
  # database_name ="GTEx v8"

  # junID <- "chr4:94334678-94342794:-" 
  # transcript_ENS = "ENST00000295256"
  # database_name ="GTEx v8"
  # junType <- 'Novel Acceptor'
  # 'get_transcript_to_plot': 'chr4:94334678-94342794:-' '' 'ENST00000295256' 'Novel Acceptor' 'FALSE
  
  message(junType)
  message("Visualise transcript: ",junID, " - ", transcript_ENS, " - ", database_name)
  
  
  ## Get junction color depending on the type of junction selected
  if (junType == "Annotated Intron") {
    range_color = "blue"
  } else if (junType == "Novel Acceptor") {
    range_color = "#64037d"
  } else {
    range_color = "#35B779FF"
  }
  
  
  ## Get coordinates from the junction selected
  junID_coordinates <- get_genomic_coordinates(junID)
  
  
  ## Get coordinates from the transcript selected
  hg38_transcript_to_plot <- get_transcript_to_plot( junID = junID, 
                                                     transcript.id = transcript_ENS, 
                                                     jxn.type = junType,
                                                     geneName = NULL, 
                                                     multiple = F)
  
  
  
  #############################
  ## PLOT
  #############################
  
  if ( is.null(hg38_transcript_to_plot) ) {
    
    visualise_empty_plot(text = paste0("Sorry, the selected ", str_to_lower(junType), " is not part of the MANE\ntranscript structure '",transcript_ENS,"'.")) %>% return()
    
    
  } else if ( nrow(hg38_transcript_to_plot$cds) == 0 | nrow(hg38_transcript_to_plot$utr) == 0 ) {
    
    visualise_empty_plot(text = paste0("Sorry, it has not been possible to plot the selected junction.\nIt is part from a transcript structure that does not have CDS or UTR regions defined.\nTranscript biotype (Ensembl v111): '", 
                                       str_replace(string = hg38_transcript_to_plot$exons$transcript_biotype %>% unique,
                                                   pattern = "_", replacement = " "), "'.")) %>%
      return()
    
  } else {
  
    hg38_transcript_to_plot$exons <- hg38_transcript_to_plot$exons %>% arrange(start)

    exons_to_zoom <- get_exons_to_zoom(jun.type = junType, 
                                       intron.to.zoom = junID_coordinates,
                                       transcript.to.plot = hg38_transcript_to_plot)
    
    
    ## PLOT
    
    if ( hg38_transcript_to_plot$exons %>% nrow() > 0 ) {
      
      missplicing_plot <- hg38_transcript_to_plot$exons %>%
        ggplot(aes(
          xstart = start,
          xend = end,
          y = hg38_transcript_to_plot$exons$transcript_id %>% unique()
        )) +
        ggtranscript::geom_range(
          data = hg38_transcript_to_plot$cds,
          fill = "#999999"
        ) +
        ggtranscript::geom_range(
          data = hg38_transcript_to_plot$utr,
          fill = "#333333",
          height = 0.25
        ) +
        ggtranscript::geom_intron(
          data = ggtranscript::to_intron(hg38_transcript_to_plot$exons, "transcript_id"),
          aes(strand = hg38_transcript_to_plot$exons$strand %>% unique)
        ) + 
        ggtranscript::geom_junction(
          data = junID_coordinates,
          color = range_color,
          ncp = 100, 
          junction.y.max = 0.5 
        ) +
        theme_light(base_size = 16) +
        ggforce::facet_zoom(xlim = c((min(exons_to_zoom$start)-100):(max(exons_to_zoom$end)+100))) +
        theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
              legend.position = "top",
              legend.text = element_text(size = "11"),
              legend.title = element_text(size = "11")) +
        xlab(paste0("Genomic position (chr", hg38_transcript_to_plot$seqnames %>% unique ,")")) + 
        ylab("") +
        guides(color = guide_legend(title = "Junction type: ")) +
        labs(title = paste0("Coordinates: '", junID_coordinates %>% pull(ID) %>% unique, "'"))
      
      missplicing_plot %>%
        return()
      
      
    } else {
      visualise_empty_plot() %>% return()
    }
  }
}



visualise_multiple_transcripts <- function(junID,
                                           transcript_ENS,
                                           geneName,
                                           database_name,
                                           junType) {
  
  
  # junID = "chr4:99336742-99342794:-"
  # junType <- "Novel Acceptor"
  # transcript_ENS = "ENST00000295256"
  # database_name = "all"
  # geneName <- "ADH1C"
  
  # junID = "chr10:87880439-87925512:+"
  # junType <- "Annotated Intron"
  # transcript_ENS = "ENST00000371953"
  # database_name = "TCGA"
  # geneName <- "PTEN"
  
  # junID <- "chr1:154597977-154598401:-"
  # transcript_ENS <- "ENST00000368474"
  # geneName <- "ADAR"
  # database_name = "TCGA"
  # junType <- "Annotated Intron"
  
  
  # junID <- "chr19:11458612-11463495:-"
  # transcript_ENS <- "ENST00000359227"
  # database_name = "TCGA"
  # junType <- "Novel Donor"
  # geneName <- "ELAVL3"
  
  
  
  ## Get junction color depending on the type of junction selected
  if (junType == "Annotated Intron") {
    range_color = "blue"
  } else if (junType == "Novel Acceptor") {
    range_color = "#64037d"
  } else {
    range_color = "#35B779FF"
  }
  
  
  ## Get coordinates from the transcript selected
  hg38_transcripts_to_plot <- get_transcript_to_plot( junID, 
                                                      geneName = geneName, 
                                                      transcript.id =  transcript_ENS, 
                                                      jxn.type = junType,
                                                      multiple = T) 
  
  
  
  
  #############################
  ## PLOT
  #############################
  
  if ( is.null(hg38_transcripts_to_plot) | nrow(hg38_transcripts_to_plot$cds) == 0 | nrow(hg38_transcripts_to_plot$utr) == 0 ) {
    
    visualise_empty_plot(text = paste0("Sorry, it has not been possible to plot the selected ", str_to_lower(junType), ".\nIt is part of a transcript structure that does not have CDS or UTR regions defined.\nTranscript biotype (Ensembl v111): '", 
                                       str_replace_all(string = hg38_transcripts_to_plot$exons$transcript_biotype %>% unique, pattern = "_", replacement = " "), "'.")) %>% return()
    
  } else {
    
    transcripts_to_plot <- hg38_transcripts_to_plot$cds %>% 
      filter(is.na(tag)) %>%
      pull(transcript_id) %>% 
      unique %>% 
      sort
    
    exons_to_plot <- hg38_transcripts_to_plot$exons %>%
      filter(transcript_id %in% transcripts_to_plot) %>%
      arrange(desc(transcript_id))
    
    ## Get coordinates from the junction selected and add transcript ID info
    junction_to_plot <- cbind(data.frame(transcript_id = exons_to_plot$transcript_id %>% unique),
                              get_genomic_coordinates(junID))
    
    ## PLOT
    if ( transcripts_to_plot %>% length() > 0 ) {
      
      missplicing_plot <- exons_to_plot %>%
        ggplot(aes(
          xstart = start,
          xend = end,
          y = transcript_id 
        )) +
        geom_range(
          data = hg38_transcripts_to_plot$cds,
          fill = "#999999"
        ) +
        geom_range(
          data = hg38_transcripts_to_plot$utr,
          fill = "#333333",
          height = 0.25
        ) +
        ggtranscript::geom_intron(
          data = ggtranscript::to_intron(exons_to_plot, "transcript_id"),
          aes(strand = ggtranscript::to_intron(exons_to_plot, "transcript_id")$strand %>% unique)
        ) + 
        ggtranscript::geom_junction(
          data = junction_to_plot,
          aes(y = transcript_id),
          color = range_color,
          ncp = 100, 
          junction.y.max = 0.5 
        ) +
        theme_light(base_size = 16) +
        theme(#axis.text.y = element_text(angle = 45, hjust = 0.5),
              legend.position = "top",
              legend.text = element_text(size = "11"),
              legend.title = element_text(size = "11")) +
        xlab(paste0(exons_to_plot$seqnames %>% unique ," (hg38)")) + 
        ylab("") +
        guides(color = guide_legend(title = "Junction type: ")) +
        labs(subtitle = paste0("Junction coordinates: '", junction_to_plot %>% pull(ID) %>% unique, "'"))
      
      missplicing_plot %>%
        return()
      
      
    } else {
      visualise_empty_plot() %>% return()
    }
  }
}

visualise_CLIP <- function(junID, 
                           transcript_ENS = NULL, 
                           geneName, 
                           junType) {

  
  # junID <- 'chr4:94334678-94342794:-'
  # geneName <- 'HPGDS'
  # junType = 'Novel Acceptor'
  # transcript_ENS = "ENST00000295256"
  
  # junID <- 'chr6:89603110-89605849:-'
  # geneName <- 'LYRM2'
  # junType = 'Novel Donor'
  # transcript_ENS = "ENST00000369393"
  
  # junID <- 'chr4:89828185-89828473:-'
  # junID <- 'chr4:89726661-89745540:-'
  # geneName <- 'SNCA' 
  # junType = 'Novel Donor'
  # transcript_ENS = NULL

  # junID <- "chr1:154597977-154598075:-"
  # geneName <- "ADAR"
  # transcript_ENS = "ENST00000368474"
  # junType ="Novel Donor"
  
  
  # junID <- "chr17:45894687-45943908:+"
  # geneName <- "MAPT"
  # transcript_ENS = "ENST00000262410"
  # junType ="Novel Acceptor"
  
  # junID <- "chr19:4501318-4501913:+"
  # geneName <- "HDGFL2"
  # transcript_ENS = "ENST00000616600"
  # junType ="Annotated Intron"
  
  # junID <- "chr19:17621871-17623117:-"
  # geneName <- "UNC13A"
  
  # junID <- "chr1:154597977-154598401:-"
  
  # junID <- "chr1:154597977-154598075:-"
  # geneName <- "ADAR"
  
  # junID <- "chr10:87862643-87864453:+"
  # transcript_ENS <- "ENST00000688308"
  # geneName <- "PTEN"
  # junType <- "Novel Donor"
  
  
  # junID <- "chr4:89835693-89836742:-"
  # geneName <- "SNCA"
  # transcript_ENS = ""
  # junType ="Annotated Intron"
  
  # junID <- "chr4:99930622-99939216:-"
  # geneName <- "DNAJB14"
  # transcript_ENS = ""
  # junType ="Novel Donor"
  
  # junID <- "chr4:94260634-94302145:-"
  # geneName <- "HPGDS"
  # transcript_ENS = "ENST00000295256"
  # junType ="Novel Acceptor"
  
  ## Get junction color depending on the type of junction selected
  if (junType == "Annotated Intron") {
    range_color = "blue"
  } else if (junType == "Novel Acceptor") {
    range_color = "#64037d"
  } else {
    range_color = "#288a5b"
  }
  
  print(paste0("CLIP: '", junID, "' '", transcript_ENS, "' - '", (geneName )[1], "' - '", junType, "'"))
  
  junID_coordinates <- get_genomic_coordinates(junID)
  
  
  ## Connect to the CLIP data database and retrieve data from the current gene
  database_path <- file.path("./database/clip_data.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  query <- paste0("SELECT * FROM '", (geneName)[1], "'")
  db_clip_data_gr <- DBI::dbGetQuery(con, query) %>%
    arrange(desc(totalClipExpNum)) %>%
    #filter(totalClipExpNum >= 10) %>%
    GenomicRanges::GRanges()
  DBI::dbDisconnect(conn = con)
  
  
  ## Get overlaps
  encori_overlaps <- GenomicRanges::findOverlaps(query = junID_coordinates %>% GenomicRanges::GRanges(),
                                                 subject = db_clip_data_gr,
                                                 maxgap = 100,
                                                 type = "any")
  
  if ( encori_overlaps %>% length() > 0 ) {
    
    ## GET MANE TRANSCRIPT INFO
    #hg38_transcript_to_plot <- get_transcript_to_plot( junID = junID, transcript.id = transcript_ENS, jxn.type = junType)
    
    ## Get coordinates from the selected transcript
    hg38_transcripts_to_plot <- get_transcript_to_plot( junID = junID, 
                                                        transcript.id = transcript_ENS, 
                                                        jxn.type = junType, 
                                                        geneName = geneName,
                                                        multiple = T ) 
    
    
    
    ## PLOT USING TRANSCRIPT INFO
    if ( is.null(hg38_transcripts_to_plot) ) {
      
      visualise_empty_plot(text = paste0("Sorry, it has not been possible to plot the CLIP-seq information for the selected ", str_to_lower(junType), 
                                         ".\nIt is not part of the MANE transcript structure.")) %>% return()
      
      
    } else if ( nrow(hg38_transcripts_to_plot$cds) == 0 | nrow(hg38_transcripts_to_plot$utr) == 0 ) {
      
      visualise_empty_plot(text = paste0("Sorry, it has not been possible to plot the CLIP-seq information for the selected junction.\nIt is part from a transcript structure that does not have CDS or UTR regions defined.\nTranscript biotype (Ensembl v111): '", 
                                         str_replace(string = hg38_transcript_to_plot$exons$transcript_biotype %>% unique,
                                                     pattern = "_", replacement = " "), "'.")) %>%
        return()
      
    } else {
      
      
      if ( any(hg38_transcripts_to_plot$transcript$transcript_id %>% unique %in% transcript_ENS) ) {
        hg38_transcript_ID_to_plot <- transcript_ENS
      } else {
        hg38_transcript_ID_to_plot <- (hg38_transcripts_to_plot$utr$transcript_id %>% unique)[1]
      }
      
      hg38_transcript_to_plot <- list(transcript = hg38_transcripts_to_plot$transcript %>% filter(transcript_id == hg38_transcript_ID_to_plot),
                                      exons = hg38_transcripts_to_plot$exons %>% filter(transcript_id == hg38_transcript_ID_to_plot),
                                      cds = hg38_transcripts_to_plot$cds %>% filter(transcript_id == hg38_transcript_ID_to_plot), 
                                      utr = hg38_transcripts_to_plot$utr %>% filter(transcript_id == hg38_transcript_ID_to_plot))
    
      
      ## GET UNIQUE TRANSCRIPT SITES
      RBP_CLIP_sites_to_plot <- db_clip_data_gr[S4Vectors::subjectHits(encori_overlaps), ] %>% 
        as_tibble() %>%
        dplyr::group_by(seqnames, start, end) %>%
        distinct(RBP, .keep_all = T) %>%
        ungroup() %>% 
        as_tibble(rownames = "index") %>% 
        mutate(index = index %>% as.integer())
      
      if ( RBP_CLIP_sites_to_plot %>% nrow > 50 ) {
        RBP_CLIP_sites_to_plot <- RBP_CLIP_sites_to_plot[1:50,]
      }
      
      
      
      hg38_transcript_to_plot$exons <- hg38_transcript_to_plot$exons %>% arrange(start)
      exons_to_zoom <- get_exons_to_zoom(jun.type = junType, 
                                         intron.to.zoom = junID_coordinates,
                                         transcript.to.plot = hg38_transcript_to_plot)
      exons_to_zoom <- exons_to_zoom %>% mutate(start = start - 100, end = end + 100)
      
      
   
      
      ## PLOT
      plot_double_annotated <- hg38_transcript_to_plot$exon %>%
        ggplot(aes(
          xstart = start,
          xend = end,
          y = hg38_transcript_to_plot$exon$transcript_id %>% unique() 
        )) +
        ggtranscript::geom_intron(
          data = to_intron(hg38_transcript_to_plot$exon, "transcript_id"),
          aes(strand = hg38_transcript_to_plot$exon$strand %>% unique)
        )  +
        geom_range(
          data = hg38_transcript_to_plot$cds,
          fill = "#cccccc"
        ) +
        geom_range(
            data = hg38_transcript_to_plot$utr,
            fill = "#333333",
            height = 0.25
            )     +
        geom_junction(
          data = junID_coordinates %>% as_tibble(),
          colour = range_color,
          junction.y.max = 0.5 ,
          ncp = 200
        )+
        ggforce::facet_zoom(xlim = c((min(exons_to_zoom$start)-100):(max(exons_to_zoom$end)+100))) +
        theme_light(base_size = 16) +
        theme(axis.text.y = element_text(angle = 90, hjust = 0.5,color = "white"),
              legend.title = element_text(size = "12"),
              legend.text  = element_text(size = "11"),
              plot.caption = element_text(size = "12"),
              plot.subtitle = element_text(size = "15"),
              legend.position = "right")+
        ylab("") +
        ylab(hg38_transcript_to_plot$exon$transcript_id %>% unique) +
        theme(axis.text.x = element_blank(),
              plot.margin = margin(t = 2, r = 2, b = -2, l = 2,unit = "pt") )

     
      clip_plot <- ggplot(data = RBP_CLIP_sites_to_plot) +
        geom_rect(
                  mapping = aes(xmin = start, xmax = end,
                                ymin = 2*index, ymax = 2*(index+1), 
                                fill = RBP)) +
        geom_label(aes(x=start+(end-start)/2, 
                      y=(2*index)+((2*(index+1))-(2*index))/2,
                      label=factor(RBP),
                      fill = factor(RBP)), 
                   colour = "white", fontface = "bold", size=3) +
                  #aes(x = start-50, y=2*(index+1), label=RBP), size=3) +
        coord_cartesian(xlim=c(min(exons_to_zoom$start)-100,
               max(exons_to_zoom$end)+100))  +
        guides(fill = guide_legend(nrow = 2)) +
        theme_light(base_size = 16) +
        theme(axis.text.y = element_text(angle = 90, hjust = 0.5,color = "white"),
              legend.title = element_text(size = "12"),
              legend.text  = element_text(size = "11"),
              plot.caption = element_text(size = "12"),
              plot.subtitle = element_text(size = "15"))+
        theme(legend.position = "none",
              plot.margin = margin(t = 2, r = 2, b = 1, l = 2,unit = "pt")) +
        ylab(hg38_transcript_to_plot$exon$transcript_id %>% unique)+
        xlab(paste0(RBP_CLIP_sites_to_plot$seqnames %>% unique() ," (hg38)")) +
        labs(caption = "* Regulatory data sourced from: (Li JH, et al. Nucleic Acids Res. 2014 Jan;42:D92-7)")# and\n(Zhou KR, Liu S, Li B, Liu SR, Zheng WJ, Cai L, et al. An encyclopedia of RNA interactomes in ENCORI).")
      
      
      ggpubr::ggarrange( plot_double_annotated,
                         clip_plot,
                         ncol = 1,
                         nrow = 2,
                         heights = c(1,2)) +
        ylab(hg38_transcript_to_plot$exon$transcript_id %>% unique )  %>%
        return()
      
    } 
    
  } else {
    visualise_empty_plot(text = paste0("Sorry, there is no CLIP-seq data available for gene  '", geneName,"'\nsupported by at least 10 CLIP experiments.\nClick on 'Download all regulatory data' to download all CLIP-seq data\nfor gene '", geneName,"'.")) %>% return()
   
  }
  
}


visualise_clinvar <- function(junID, 
                              #clinvar_locus,
                              
                              #CLNSIG_list,
                              #CLNVC_list,
                              #MC_list,
                              
                              database_name, 
                              junType, 
                              geneName) {
  
  # junID = "chr1:155235308-155235680:-"
  # database_name = "GTEx v8" 
  # junType = "Novel Acceptor" 
  # geneName = "GBA1" 
 
  # junID = "chr1:154585345-154585894:-"
  # database_name = "GTEx v8" 
  # junType = "Novel Donor" 
  # geneName = "ADAR"
  
  # junID = "chr10:87862643-87925512:+"
  # geneName = "PTEN"
  # junType = "Novel Donor" 
  # database_name = "GTEx v8"
  
  # junID = 'chr4:89822344-89828142:-' 
  # database_name = 'TCGA' 
  # junType = 'Novel Acceptor' 
  # geneName = 'SNCA'
  
  message("CLINVAR visualisation: '", junID, "' '", database_name, "' '", junType, "' '", geneName, "'")
  
  ## Get junction color depending on the type of junction selected
  if (junType == "Annotated Intron") {
    range_color = "blue"
  } else if (junType == "Novel Acceptor") {
    range_color = "#64037d"
  } else {
    range_color = "#288a5b"
  }
  
  ## Get junction coordinates
  junID_coordinates <- get_genomic_coordinates(junID)
    
  ## Get clinvar data
  clinvar_data <- get_database_clinvar(junID_coordinates) 
  
  clinvar_data <- clinvar_data %>% 
    mutate(strand = junID_coordinates$strand %>% as.character())
  
  ## Get the transcript
  #hg38_transcript_to_plot <- get_transcript_to_plot(junID = junID, geneName = geneName, jxn.type = junType)
  
  ## Get coordinates from the transcript selected
  hg38_transcript_to_plot <- get_transcript_to_plot( junID = junID, 
                                                      transcript.id = NULL, 
                                                      jxn.type = junType, 
                                                      geneName = geneName,
                                                      multiple = F) 

  ## Plot
  ## PLOT USING TRANSCRIPT INFO
  if ( is.null(hg38_transcript_to_plot) ) {
    
    visualise_empty_plot(text = paste0("Sorry, it has not been possible to plot the CLINVAR information for the selected ", str_to_lower(junType), 
                                       ".\nIt is not part of a MANE transcript structure.")) %>% return()
    
    
  } else if ( nrow(hg38_transcript_to_plot$cds) == 0 | nrow(hg38_transcript_to_plot$utr) == 0 ) {
    
    visualise_empty_plot(text = paste0("Sorry, it has not been possible to plot the CLINVAR information for the selected junction.\nIt is part from a transcript structure that does not have CDS or UTR regions defined.\nTranscript biotype (Ensembl v111): '", 
                                       str_replace_all(string = hg38_transcript_to_plot$exons$transcript_biotype %>% unique,
                                                       pattern = "_", replacement = " "), "'.")) %>%
      return()
    
  } else {
    
    
    ## GET THE EXONS TO ZOOM IN
    hg38_transcript_to_plot$exons <- hg38_transcript_to_plot$exons %>% arrange(start)
    
    exons_to_zoom <- get_exons_to_zoom(jun.type = junType, 
                                       intron.to.zoom = junID_coordinates,
                                       transcript.to.plot = hg38_transcript_to_plot)
    
    
    ## PLOT
    plot_double_annotated <- hg38_transcript_to_plot$exons %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = hg38_transcript_to_plot$exons$transcript_id %>% unique() 
      )) +
      ggtranscript::geom_intron(
        data = to_intron(hg38_transcript_to_plot$exons, "transcript_id"),
        aes(strand = hg38_transcript_to_plot$exons$strand %>% unique)
      ) +   
      geom_range(
        data = hg38_transcript_to_plot$cds,
        fill = "#cccccc"
      ) +
      geom_range(
        data = hg38_transcript_to_plot$utr,
        fill = "#333333",
        height = 0.25
      ) +
      
      ggrepel::geom_label_repel(
        data = clinvar_data,
        aes(x = start, 
            fontface = "bold",
            label = paste0(ID),
            color = CLNSIG),
        #nudge_y = 0.5,
        #nudge_x = 3,
        #direction = c("y"),
        max.overlaps = 100
      ) +
      
      ggtranscript::geom_junction(
        data = junID_coordinates,
        ncp = 100, 
        colour = range_color,
        junction.y.max = 0.5 
      ) +
      
      ggforce::facet_zoom(xlim = c((min(exons_to_zoom$start)-100):(max(exons_to_zoom$end)+100))) +
      theme_light(base_size = 16) +
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
            legend.title = element_text(size = "10"),  
            legend.text  = element_text(size = "10"), 
            plot.caption = element_text(size = "12"),
            plot.subtitle = element_text(size = "15"),
            legend.position = "top") +
      #guides(fill = guide_legend(element_blank())) +
      xlab(paste0("Genomic position (", hg38_transcript_to_plot$exons$seqnames %>% unique() ,")")) + 
      ylab(paste0(geneName)) +
      labs(#subtitle = "Displaying regulatory data from 'HepG2' and 'K562' cell lines:",
           caption = "*ClinVar data sourced from: (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/).")
    
    plot_double_annotated
    
    ## Only annotated the zoomed panel
    plot_zoomed_annotated <- ggplot_build(plot_double_annotated)
    plot_zoomed_annotated$data[[4]][which( plot_zoomed_annotated$data[[4]][,"PANEL"] == 1 ),"alpha"] <- 0
    plot_zoomed_annotated <- ggplot_gtable(plot_zoomed_annotated)
    plot(plot_zoomed_annotated) %>%
      return()
    
  } 
}


