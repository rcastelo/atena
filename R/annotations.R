#' Get RepeatMasker UCSC annotations
#'
#' The \code{annotaTEs()} function fetches RepeatMasker UCSC transposable
#' element (TE) annotations using 
#' \link[AnnotationHub]{AnnotationHub} and parses them.
#' 
#' @param genome The genome version of the desired RepeatMasker annotations
#'        (e.g. "hg38").
#'
#' @param parsefun A \code{function} to parse the annotations:
#' \itemize{
#'   \item Function \code{rmskidentity} returns RepeatMasker annotations as
#'         present in \link[AnnotationHub]{AnnotationHub}, without
#'         processing them.
#'   \item Function \code{rmskbasicparser} parses annotations by removing
#'         low complexity regions, simple repeats, satellites, rRNA, scRNA,
#'         snRNA, srpRNA and tRNA. Also removes TEs
#'         with a strand different than "+" or "-". Modifies "repFamily" and
#'         "repClass" columns when a "?" is present or when they are defined
#'         as "Unknown" or "Other". Finally, assigns a unique id to each TE
#'         instance by adding the suffix "_dup" plus a number at the end of
#'         the "repName".
#'   \item User-defined function. Input and output should be \code{GRanges}
#'         objects.
#' }
#'        
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object with
#'         transposable element annotations.
#'
#' @details
#' Given a specific genome version, the \code{annotaTEs()} function fetches
#' RepeatMasker annotations from UCSC Genome Browser using the 
#' \link[AnnotationHub]{AnnotationHub} package. Since RepeatMasker not only
#' provides TE annotations but also low complexity DNA sequences and other
#' types of repeats, a specific \code{parsefun} can be set to parse these
#' annotations (e.g. \code{rmskbasicparser} or a user-defined function). If no
#' parsing is required, \code{parsefun} can be set to \code{rmskidentity}.
#'
#' @seealso
#' \code{\link[AnnotationHub]{AnnotationHub}}
#' 
#'
#' @examples
#' rmsk_gr <- annotaTEs(genome = "hg38", parsefun = rmskidentity)
#' 
#' 
#' @aliases annotaTEs
#' @rdname annotaTEs
#' @name annotaTEs
#' @importFrom AnnotationHub AnnotationHub query
#' @export
annotaTEs <- function(genome="hg38", parsefun=rmskidentity) {
  ah <- AnnotationHub()
  qah <- query(ah, c(genome, "RepeatMasker", "UCSC"))
  if (length(qah) == 0)
    stop(sprintf("UCSC RepeatMasker tracks for genome %s not found", genome))
  else if (length(qah) > 1)
    stop(sprintf("more than one UCSC RepeatMasker track for genome %s found", genome))
  
  id <- names(qah)
  gr <- ah[[id]]
  parsefun(gr)
}

#' Parser of RepeatMasker annotations
#' @param gr A \link[GenomicRanges:GRanges-class]{GRanges} object with
#'           RepeatMasker annotations from \link[AnnotationHub]{AnnotationHub}
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object.
#'         
#' @details 
#' Parses annotations by removing low complexity regions, simple repeats,
#' satellites, rRNA, scRNA, snRNA, srpRNA and tRNA. Also removes TEs with a
#' strand different than "+" or "-". Modifies "repFamily" and
#' "repClass" columns when a "?" is present or when they are defined
#' as "Unknown" or "Other". Finally, assigns a unique id to each TE instance
#' by adding the suffix "_dup" plus a number at the end of the "repName".
#'
#' @examples
#' rmsk_gr <- annotaTEs(genome = "dm6", parsefun = rmskbasicparser)
#' 
#' @aliases rmskbasicparser
#' @rdname rmskbasicparser
#' @name rmskbasicparser
#' @importFrom GenomicRanges strand
#' @export
rmskbasicparser <- function(gr) {
  to_discard <- c("Low_complexity", "Simple_repeat", "Satellite", "rRNA",
                  "scRNA", "snRNA", "srpRNA", "tRNA")
  gr <- gr[!(gr$repClass %in% to_discard),]
  gr <- gr[strand(gr) == "+" | strand(gr) == "-",]
  # removing '?' from name, class or family
  gr$repName <- gsub(pattern ="?", x =gr$repName, replacement="", fixed =TRUE)
  gr$repClass <- gsub(pattern ="?", x=gr$repClass, replacement="", fixed=TRUE)
  gr$repFamily <- gsub(pattern ="?", x =gr$repFamily, replacement="", 
                       fixed =TRUE)
  gr <- unlist(split(gr, f = gr$repName))
  dupl_te <- duplicated(gr$repName) | duplicated(gr$repName, fromLast = TRUE)
  ndup <- lengths(split(gr, f = gr$repName))
  dup <- unlist(lapply(ndup, FUN=function(x) seq_len(x)))
  repName <- paste(gr$repName, "_dup", dup, sep ="")
  repName[!dupl_te] <- gr$repName[!dupl_te]
  names(gr) <- repName
  gr <- sort(gr)

  i <- gr$repFamily %in% c("Unknown", "Other")
  if (any(i)) {
    gr$repFamily[i] <- gr$repName[i]
  }
  gr
}

#' Identity function for parsefun
#' @param gr A \link[GenomicRanges:GRanges-class]{GRanges} object.
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object.
#'
#' @details 
#' Identity function: returns the \link[GenomicRanges:GRanges-class]{GRanges}
#' object without any modification.
#' 
#' @examples
#' rmsk_gr <- annotaTEs(genome = "dm6", parsefun = rmskidentity)
#'         
#' @aliases rmskidentity
#' @rdname rmskidentity
#' @name rmskidentity
#' @export
rmskidentity <- function(gr) {
  gr
}


