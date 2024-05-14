
visualise_transcript <- function(junID,
                                 transcript_ENS,
                                 database_name,
                                 junType) {
  
  
  # junID = "chr1:154597977-154598075:-" 
  # transcript_ENS = "ENST00000647682"
  # database_name = "LRRK2 G2019S mutation"
  # junType <- "Novel Donor"
  
  # junType <- "Novel Donor"
  # junID = "chr10:87862643-87864453:+"
  # transcript_ENS = "ENST00000688308"
  # database_name = "TCGA"
  
  message(junType)
  message(junID, " - ", transcript_ENS, " - ", database_name)
  
  
  ## Get junction color depending on the type of junction selected
  if (junType == "Annotated Intron") {
    range_color = "blue"
  } else if (junType == "Novel Acceptor") {
    range_color = "#64037d"
  } else {
    range_color = "#35B779FF"
  }
  
  
  ## Get coordinates from the junction selected
  df_intron_gr <- get_genomic_coordinates(junID) %>% GenomicRanges::GRanges()
  
  
  ## Get coordinates from the transcript selected
  hg38_transcript_to_plot <- get_transcript_to_plot( junID = junID, transcript_id = transcript_ENS ) 
  
  
  #############################
  ## PLOT
  #############################
  
  if ( is.null(hg38_transcript_to_plot) | nrow(hg38_transcript_to_plot$cds) == 0 | nrow(hg38_transcript_to_plot$utr) == 0 ) {
    
    ggplot() +
      theme_void(base_size = 14) +
      geom_text(aes(0,0,
                    label=paste0("Sorry, it has not been possible to plot the selected junction.\nThe transcript might not have CDS or UTR regions defined.\nTranscript biotype (Ensembl v111): '", 
                                 str_replace(string = hg38_transcript_to_plot$exons$transcript_biotype %>% unique,
                                             pattern = "_", replacement = " "), "'."))) %>%
      return()
    
  } else {
  

    ## GET THE EXONS THAT WILL BE ZOOMED IN
    index_first_exon <- which(abs(hg38_transcript_to_plot$exons$start - (df_intron_gr %>% GenomicRanges::end())) == 
                                min(abs(hg38_transcript_to_plot$exons$start-(df_intron_gr %>% GenomicRanges::end()))))
    index_second_exon <- which(abs(hg38_transcript_to_plot$exons$end - (df_intron_gr %>% GenomicRanges::start())) == 
                                 min((abs(hg38_transcript_to_plot$exons$end - (df_intron_gr %>% GenomicRanges::start())))[-index_first_exon]))
    
    exons_to_zoom <- rbind(hg38_transcript_to_plot$exons[index_first_exon,],
                           hg38_transcript_to_plot$exons[index_second_exon,])
    
    print(hg38_transcript_to_plot$exons$strand %>% unique)
    ## PLOT
    
    if ( hg38_transcript_to_plot$exons %>% nrow() > 0 ) {
      
      missplicing_plot <- hg38_transcript_to_plot$exons %>%
        ggplot(aes(
          xstart = start,
          xend = end,
          y = hg38_transcript_to_plot$exons$transcript_id %>% unique()
        )) +
        ggtree::geom_range(
          data = hg38_transcript_to_plot$cds,
          fill = "#999999"
        ) +
        ggtree::geom_range(
          data = hg38_transcript_to_plot$utr,
          fill = "#333333",
          height = 0.25
        ) +
        ggtranscript::geom_intron(
          data = ggtranscript::to_intron(hg38_transcript_to_plot$exons, "transcript_id"),
          aes(strand = hg38_transcript_to_plot$exons$strand %>% unique)
        ) + 
        ggtranscript::geom_junction(
          data = df_intron_gr %>% as_tibble(),
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
        labs(title = paste0("Coordinates: '", df_intron_gr %>% as_tibble() %>% pull(ID) %>% unique, "'"))
      
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
  hg38_transcripts_to_plot <- get_transcript_to_plot( junID, geneName = geneName, transcript_id =  transcript_ENS, multiple = T) 
  
  
  
  
  #############################
  ## PLOT
  #############################
  
  if ( is.null(hg38_transcripts_to_plot) | nrow(hg38_transcripts_to_plot$cds) == 0 | nrow(hg38_transcripts_to_plot$utr) == 0 ) {
    
    ggplot() +
      theme_void(base_size = 14) +
      geom_text(aes(0,0,
                    label=paste0("Sorry, it has not been possible to plot the selected junction.\nThe transcript might not have CDS or UTR regions defined.\nTranscript biotype (Ensembl v111): '", 
                                 str_replace(string = hg38_transcript_to_plot$exons$transcript_biotype %>% unique,
                                             pattern = "_", replacement = " "), "'."))) %>%
      return()
    
  } else {
    
    transcripts_to_plot <- hg38_transcripts_to_plot$cds %>% pull(transcript_id) %>% unique %>% sort
    
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
          aes(strand = hg38_transcripts_to_plot$strand )
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

visualise_CLIP <- function(junID, transcript_ENS, geneName, junType) {
  
  # print(junID, " - ", transcript_ENS, " - ", geneName, " - ", junType)
  # junID <- "chr19:17621871-17623117:-"
  # geneName <- "UNC13A"
  # junID <- "chr1:154597977-154598401:-"
  # junID <- "chr1:154597977-154598075:-"
  # geneName <- "ADAR"
  
  # junID <- "chr10:87862643-87864453:+"
  # transcript_ENS <- "ENST00000688308"
  # geneName <- "PTEN"
  # junType <- "Novel Donor"
  
  ## Get junction color depending on the type of junction selected
  if (junType == "Annotated Intron") {
    range_color = "blue"
  } else if (junType == "Novel Acceptor") {
    range_color = "#64037d"
  } else {
    range_color = "#288a5b"
  }
  
  print(paste0("CLIP: ", junID, " ", transcript_ENS, " - ", (geneName )[1]))
  
  junID_gr <- get_genomic_coordinates(junID) %>% 
    GenomicRanges::GRanges()
  
  ## Connect to the CLIP data database and retrieve data from the current gene
  database_path <- file.path("./database/clip_data.sqlite")
  
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  query <- paste0("SELECT * FROM '", (geneName )[1], "'")
  db_clip_data_gr <- DBI::dbGetQuery(con, query) %>%
    filter(totalClipExpNum >= 10) %>%
    GenomicRanges::GRanges()
  
  DBI::dbDisconnect(conn = con)
  
  ## Get overlaps
  encori_overlaps <- GenomicRanges::findOverlaps(query = junID_gr,
                                                 subject = db_clip_data_gr,
                                                 maxgap = 100,
                                                 type = "any")
  
  if ( encori_overlaps %>% length() > 0 ) {
    
    ## GET MANE TRANSCRIPT INFO
    hg38_transcript_to_plot <- get_transcript_to_plot( junID = junID, transcript_id = transcript_ENS)
    
    ## PLOT USING TRANSCRIPT INFO
    if ( !is.null(hg38_transcript_to_plot) ) {
      
      ## GET UNIQUE TRANSCRIPT SITES
      RBP_CLIP_sites <- db_clip_data_gr[S4Vectors::subjectHits(encori_overlaps), ] %>% 
        as_tibble() %>%
        dplyr::group_by(seqnames,start,end) %>%
        distinct(RBP, .keep_all = T) %>%
        ungroup()
      
      
      ## GET THE EXONS THAT WILL BE ZOOMED IN
      index_first_exon <- which(abs(hg38_transcript_to_plot$exons$start - (junID_gr %>% GenomicRanges::end())) == min(abs(hg38_transcript_to_plot$exons$start-(junID_gr %>% GenomicRanges::end()))))
      index_second_exon <- which(abs(hg38_transcript_to_plot$exons$end - (junID_gr %>% GenomicRanges::start())) == min((abs(hg38_transcript_to_plot$exons$end - (junID_gr %>% GenomicRanges::start())))[-index_first_exon]))

      exons_to_zoom <- rbind(hg38_transcript_to_plot$exons[index_first_exon,],
                             hg38_transcript_to_plot$exons[index_second_exon,])
      exons_to_zoom <- exons_to_zoom %>% mutate(start = start-100, end=end+100)
      
      
   
      
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
        )  +geom_range(
          data = hg38_transcript_to_plot$cds,
          fill = "#cccccc"
        ) +
        geom_range(
            data = hg38_transcript_to_plot$utr,
            fill = "#333333",
            height = 0.25
            )     +
        geom_junction(
          data = junID_gr %>% as_tibble(),
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
      #+
        #labs(subtitle = "Regulatory data suported by at least 10 CLIP experiments:")
        
        # geom_junction_label_repel(
        #   data = junID_gr %>% as_tibble(),
        #   aes(label = ID),
        #   colour = range_color,
        #   junction.y.max = 0.5,
        #   segment.color = NA
        # )  
        # #ggforce::facet_zoom(xlim = c((min(exons_to_zoom$start)-100):(max(exons_to_zoom$end)+100))) +
        # theme_light(base_size = 16) +
        # theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        #       legend.title = element_text(size = "12"),  
        #       legend.text  = element_text(size = "11"), 
        #       plot.caption = element_text(size = "12"),
        #       plot.subtitle = element_text(size = "15"),
        #       legend.position = "right") +
        # guides(fill = guide_legend(ncol = 2))# +
        #xlab(paste0(hg38_transcript_to_plot$exons$seqnames %>% unique() ," (hg38)")) + 
        #ylab(paste0(geneName)) +
        #labs(subtitle = "Displaying regulatory data suported by more than 10 CLIP experiments:",
        #     caption = "*Regulatory data sourced from: (Li JH, et al. Nucleic Acids Res. 2014 Jan;42:D92-7) and\n(Zhou KR, Liu S, Li B, Liu SR, Zheng WJ, Cai L, et al. An encyclopedia of RNA interactomes in ENCORI).")
      
      #plot_double_annotated 
      # Only annotated the zoomed panel
      #plot_zoomed_annotated <- ggplot_build(plot_double_annotated)
      #plot_zoomed_annotated$data[[4]][which( plot_zoomed_annotated$data[[4]][,"PANEL"] == 1 ),"alpha"] <- 0
      #plot_zoomed_annotated <- ggplot_gtable(plot_zoomed_annotated)
      
      
     
      RBP_CLIP_sites_to_plot <- RBP_CLIP_sites %>% 
        arrange(start)  %>% 
        as_tibble(rownames = "index") %>% 
        mutate(index = index %>% as.integer())
      
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
      
    } else {
      visualise_empty_plot()
    }
    
  } else {
    visualise_empty_plot()
  }
  
}


visualise_clinvar <- function(junID, clinvar_locus, database_name, junType, geneName) {
  
  # clinvar_locus <- "chr1:154597977"
 
  # junID = "chr1:155235308-155235680:-"
  # clinvar_locus = c("chr1:155235101,chr1:155235193,chr1:155235679,chr1:155235680,chr1:155237579,chr1:155239615")
  # database_name = "GTEx v8" 
  # junType = "Novel Acceptor" 
  # geneName = "GBA1"
  
  message(junID, " ", clinvar_locus, " ", database_name, " ", junType, " ", geneName)
  
  ## Get junction color depending on the type of junction selected
  if (junType == "Annotated Intron") {
    range_color = "blue"
  } else if (junType == "Novel Acceptor") {
    range_color = "#64037d"
  } else {
    range_color = "#288a5b"
  }
  
  ## Get junction coordinates
  junID_coordinates_gr <- get_genomic_coordinates(junID)%>% GenomicRanges::GRanges()
    
  ## Get clinvar coordinates
  clinvar_locus_coordinates <- get_genomic_coordinates(coordinates = clinvar_locus) 
  
  ## Get the transcript
  transcript_info <- get_transcript_to_plot(junID = junID, geneName = geneName)

  ## Plot
  if ( !is.null(transcript_info) ) {
    
    ## GET THE EXONS THAT WILL BE ZOOMED IN
    index_first_exon <- which(abs(transcript_info$exons$start - (junID_coordinates_gr %>% GenomicRanges::end())) == min(abs(transcript_info$exons$start-(junID_coordinates_gr %>% GenomicRanges::end()))))
    index_second_exon <- which(abs(transcript_info$exons$end - (junID_coordinates_gr %>% GenomicRanges::start())) == min((abs(transcript_info$exons$end - (junID_coordinates_gr %>% GenomicRanges::start())))[-index_first_exon]))
    
    exons_to_zoom <- rbind(transcript_info$exons[index_first_exon,],
                           transcript_info$exons[index_second_exon,])
    
    ## PLOT
    plot_double_annotated <- transcript_info$exons %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = transcript_info$exons$transcript_id %>% unique() 
      )) +
      ggtranscript::geom_intron(
        data = to_intron(transcript_info$exons, "transcript_id"),
        aes(strand = transcript_info$exons$strand %>% unique)
      ) +   
      geom_range(
        data = transcript_info$cds,
        fill = "#cccccc"
      ) +
      geom_range(
        data = transcript_info$utr,
        fill = "#333333",
        height = 0.25
      ) +
      ggrepel::geom_label_repel(
        data = clinvar_locus_coordinates,
        aes(x = start, fontface = "bold",
            label = paste0(ID),#, " \n ", seqnames, ":", start, "-", end),#),
            color = ID),
        nudge_y = 0.5,
        nudge_x = 3,
        direction = c("y")
      ) +
      
      geom_junction(
        data = junID_coordinates_gr %>% as_tibble(),
        ncp = 100, 
        colour = range_color,
        junction.y.max = 0.5 
      ) +
      # geom_junction_label_repel(
      #   data = junID_coordinates_gr %>% as_tibble(),
      #   aes(label = ID),
      #   colour = range_color,
      #   junction.y.max = 0.5,
      #   segment.color = NA
      # ) +
      
      ggforce::facet_zoom(xlim = c((min(exons_to_zoom$start)-100):(max(exons_to_zoom$end)+100))) +
      theme_light(base_size = 16) +
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
            legend.title = element_text(size = "10"),  
            legend.text  = element_text(size = "10"), 
            plot.caption = element_text(size = "12"),
            plot.subtitle = element_text(size = "15"),
            legend.position = "none") +
      #guides(fill = guide_legend(element_blank())) +
      xlab(paste0("Genomic position (", transcript_info$exons$seqnames %>% unique() ,")")) + 
      ylab(paste0(geneName)) +
      labs(#subtitle = "Displaying regulatory data from 'HepG2' and 'K562' cell lines:",
           caption = "*ClinVar data sourced from: (https://www.ncbi.nlm.nih.gov/variation/view).")
    
    plot_double_annotated
    
    ## Only annotated the zoomed panel
    plot_zoomed_annotated <- ggplot_build(plot_double_annotated)
    plot_zoomed_annotated$data[[4]][which( plot_zoomed_annotated$data[[4]][,"PANEL"] == 1 ),"alpha"] <- 0
    plot_zoomed_annotated <- ggplot_gtable(plot_zoomed_annotated)
    plot(plot_zoomed_annotated) %>%
      return()
    
  } else {
    visualise_empty_plot()
  }
}


