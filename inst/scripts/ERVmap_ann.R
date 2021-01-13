## script to import annotations from human endogenous retroviruses (HERVs)
## published in
## Tokuyama M et al. ERVmap analysis reveals genome-wide transcription of human
## endogenous retroviruses. PNAS. 2018;115(50):12565-12572. DOI:

ERVmap_ann <- rtracklayer::import.bed("ERVs/ERVmap-master/ERVmap.bed")
seqlevels(ERVmap_ann) <-
  levels(as.factor(paste("chr",
                         GenomicAlignments::seqnames(ERVmap_ann),
                         sep = '')))
saveRDS(ERVmap_ann, file=file.path("..", "extdata", "ERVmap_ann.rds"))
