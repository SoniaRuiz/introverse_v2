download_CLIP_data <- function(junID, geneName) {
  
  # junID <- "chr1:154597977-154598401:-"
  # geneName <- "ADAR"
  iCLIP_RBP_file_path <- paste0("./dependencies/ENCORI/ENCORI_hg38_RBP_all_Gene_", geneName %>% unique,".txt")
  
  iCLIP_RBP_data_gr <- readr::read_delim(file = iCLIP_RBP_file_path,
                                         delim = "\t", col_names = T, skip = 3, show_col_types = FALSE) %>%
    dplyr::rename(seqnames = chromosome, start = narrowStart, end = narrowEnd) %>%
    GenomicRanges::GRanges()
  
  junID_gr <- get_genomic_coordinates(junID) %>% 
    GenomicRanges::GRanges()
  
  encori_overlaps <- GenomicRanges::findOverlaps(query = junID_gr,
                                                 subject = iCLIP_RBP_data_gr,
                                                 maxgap = 100,
                                                 type = "any")
  
  RBP_CLIP_sites <- iCLIP_RBP_data_gr[S4Vectors::subjectHits(encori_overlaps), ] %>% 
    as_tibble() %>%
    dplyr::group_by(seqnames, start, end) %>%
    distinct(RBP, .keep_all = T) %>%
    ungroup()
  
  RBP_CLIP_sites %>%
    return()
}
