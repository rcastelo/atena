rmskidentity <- function(gr) {
  gr
}

annotaTEs <- function(genome="hg38", parsefun=rmskidentity) {
  ah <- AnnotationHub()
  qah <- query(ah, c(genome, "RepeatMasker"))
  if (length(qah) == 0)
    stop(sprintf("UCSC RepeatMasker tracks for genome %s not found", genome))
  else if (length(qah) > 1)
    stop(sprintf("more than one UCSC RepeatMasker track for genome %s found", genome))

  id <- names(qah)
  gr <- ah[[id]]
  parsefun(gr)
}
