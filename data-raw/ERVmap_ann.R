
ERVmap_ann <- rtracklayer::import.bed("ERVs/ERVmap-master/ERVmap.bed")
seqlevels(ERVmap_ann) <-
  levels(as.factor(paste("chr",
                         GenomicAlignments::seqnames(ERVmap_ann),
                         sep = '')))

usethis::use_data(ERVmap_ann, overwrite = TRUE)
