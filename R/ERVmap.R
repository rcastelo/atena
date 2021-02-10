#' Build an ERVmap parameter object
#'
#' Build an object of the class \code{ERVmapParam}
#'
#' @param bfl A \code{BamFile} or \code{BamFileList} object, or a character
#' string vector of BAM filenames.
#'
#' @param annotations A \code{GRanges} or \code{GRangesList} object. Elements
#' in this object should have names, which will be used as a grouping factor
#' for ranges forming a common locus.
#'
#' @param annotations A 'GRanges' object.
#'
#' @param singleEnd (Default FALSE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @param strandMode (Default 1) Numeric vector which can take values 0, 1 or 2.
#'   The strand mode is a per-object switch on
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   objects that controls the behavior of the strand getter. See
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   class for further detail. If \code{singleEnd = TRUE}, then use either
#'   \code{strandMode = NULL} or do not specify the \code{strandMode} parameter.
#'
#' @param ignoreStrand (Default TRUE) A logical which defines if the strand
#' should be taken into consideration when computing the overlap between reads
#' and TEs/ERVs in the annotations. When \code{ignore_strand = FALSE}, the
#' \code{\link[GenomicAlignments]{summarizeOverlaps}} function will only
#' consider those reads selected after filtering which overlap the TE or
#' ERV on the same strand. On the contrary, when \code{ignore_strand = TRUE},
#' the \code{\link[GenomicAlignments]{summarizeOverlaps}} function will count
#' any alignment which overlaps with the element in the annotations regardless
#' of the strand. For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}}.
#'
#' @param filterUniqReads (Default TRUE) Logical value indicating whether to apply
#' the alignment filters to unique reads (TRUE) or not (FALSE). These filters,
#' which are always applied to multi-mapping reads, are optional for unique
#' reads. If TRUE, the unique reads not passing one or more filters from the
#' ERVmap pipeline, except for the "AS - XS >= 5" filter, will be discarded to
#' compute TEs expression.
#'
#' @param fragments (Default TRUE) A logical; applied to paired-end data only.
#' When \code{fragments=TRUE} (default), the read-counting method in the
#' original ERVmap algorithm will be applied, by which each mate of a paired-end
#' read is counted once, and therefore two mates mapping to the same element
#' result in adding up a count value of two. When \code{fragments=FALSE}, if the
#' two mates of a paired-end read map to the same element, they are counted as
#' a single hit and singletons, reads with unmapped pairs and other fragments
#' are not counted.
#'
#' @return A \linkS4class{ERVmapParam} object.
#'
#' @examples
#' bamfiles <- list.files(system.file("extdata", package="atena"),
#'                        pattern="*.bam", full.names=TRUE)
#' annot <- ERVmap_ann()
#' empar <- ERVmapParam(bamfiles, annot, singleEnd = TRUE)
#' empar
#'
#' @references
#' Tokuyama M et al. ERVmap analysis reveals genome-wide transcription of human
#' endogenous retroviruses. PNAS. 2018;115(50):12565-12572. DOI:
#' \url{https://doi.org/10.1073/pnas.1814589115}
#'
#' @importFrom methods is new
#' @export
ERVmapParam <- function(bfl, annotations,
                        singleEnd=FALSE,
                        ignoreStrand=TRUE,
                        strandMode=1L,
                        filterUniqReads=FALSE,
                        fragments=TRUE) {
  if (missing(bfl) || !class(bfl) %in% c("character", "BamFileList"))
    stop("argument 'bfl' should be either a string character vector of BAM file names or a 'BamFileList' object")

  if (is.character(bfl)) {
    mask <- sapply(bfl, file.exists)
    if (any(!mask))
      stop(sprintf("The following input BAM files cannot be found:\n%s",
                   paste(paste("  ", bfl), collapse="\n")))
  }
  if (!is(bfl, "BamFileList"))
    bfl <- BamFileList(bfl)

  annotationsobjname <- deparse(substitute(annotations))
  env <- parent.frame()
  if (!exists(annotationsobjname))
    stop(sprintf("input annotation object '%s' is not defined.", annotationsobjname))

  if (!is(annotations, "GRanges") && !is(annotations, "GRangesList"))
    stop(sprintf("annotations object '%s' should be either a 'GRanges' or a 'GRangesList' object.",
                 annotationsobjname))

  if (is.null(names(annotations)))
    stop(sprintf("the annotations object '%s' has no names.", annotationsobjname))

  if (is(annotations, "GRangesList"))
    annotations <- unlist(annotations)

  new("ERVmapParam", bfl=bfl, annotations=annotations, singleEnd=singleEnd,
      ignoreStrand=ignoreStrand, strandMode=as.integer(strandMode),
      filterUniqReads=filterUniqReads, fragments=fragments)
}

#' @param object A \linkS4class{ERVmapParam} object.
#'
#' @importFrom GenomeInfoDb seqlevels
#' @export
#' @aliases show,ERVmapParam-method
#' @rdname ERVmapParam-class
setMethod("show", "ERVmapParam",
          function(object) {
            cat(class(object), "object\n")
            cat(sprintf("# BAM files (%d): %s\n", length(object@bfl),
                        .pprintnames(names(object@bfl))))
            cat(sprintf("# annotations (%d): %s\n", length(object@annotations),
                        ifelse(is.null(names(object@annotations)),
                               paste("on", .pprintnames(seqlevels(object@annotations))),
                               .pprintnames(names(object@annotations)))))
            cat(sprintf("# %s, %s",
                        ifelse(object@singleEnd, "single-end", "paired-end"),
                        ifelse(object@ignoreStrand, "unstranded", "stranded")))
            if (!object@ignoreStrand)
              cat(sprintf(" (strandMode=%d)", object@strandMode))
            cat(sprintf(", %s",
                        ifelse(object@filterUniqReads, "unique-read filtering",
                               "unfiltered unique reads")))
            if (!object@singleEnd)
              cat(sprintf(", %s",
                          ifelse(object@fragments, "counting each paired-end mate",
                                 "counting both paired-end mates")))
            cat("\n")
          })

#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame
#' @export
#' @aliases qtex
#' @aliases qtex,ERVmapParam-method
#' @rdname qtex
setMethod("qtex", "ERVmapParam",
          function(x, phenodata=NULL, BPPARAM=SerialParam(progressbar=TRUE)) {
            if (!is.null(phenodata)) {
              if (nrow(phenodata) != length(x@bfl))
                stop("number of rows in 'phenodata' is different than the number of input BAM files in the input parameter object 'x'.")
              if (is.null(rownames(phenodata)))
                stop("'phenodata' has no row names.")
            }

            if (x@singleEnd)
              cnt <- bplapply(x@bfl, .qtex_ervmap_singleend, ervpar=x, BPPARAM=BPPARAM)
              
            else
              cnt <- bplapply(x@bfl, .qtex_ervmap_pairedend, ervpar=x, BPPARAM=BPPARAM)

            cntmat <- do.call("cbind", cnt)
            colnames(cntmat) <- gsub(".bam$", "", names(x@bfl))
            colData <- DataFrame(row.names=colnames(cntmat))
            if (!is.null(phenodata)) {
              colData <- phenodata
              colnames(cnt) <- rownames(colData)
            }

            SummarizedExperiment(assays=list(counts=cntmat),
                                 rowRanges=x@annotations,
                                 colData=colData)
          })


#' @importFrom Rsamtools scanBamFlag ScanBamParam yieldSize asMates
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom SummarizedExperiment SummarizedExperiment

.qtex_ervmap_singleend <- function(bf, ervpar) {

  tags_df <- .get_tags_in_BAM_singleend(bf)
  yieldSize(bf) <- 100000
  sbflags <- scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE,
                         isNotPassingQualityControls = FALSE,
                         isSupplementaryAlignment = FALSE)
  
  if (!tags_df$NH & !tags_df$XS) {
    if (!ervpar@filterUniqReads) {
      stop("Error in 'filterUniqReads = FALSE': Neither the NH tag nor the XS tag (from BWA) are provided in the BAM files. Unique reads cannot be differentiated from multi-mapping reads, therefore unique reads must be also filtered.")
    }
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"))
  } else if (tags_df$NH) {
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"),
                          tagFilter = list(NH = c(2:10000)))
  } else if (tags_df$XS) {
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"),
                          tagFilter = list(XS = c(1:10000)))
  }
  if (!tags_df$NM & !tags_df$nM) {
    stop("Neither the NM tag nor the nM tag are available in the BAM file. At least one of them is needed to apply the 2nd filter.")
  }
  
  cnt <- NULL
  open(bf)
  while (length(r <- readGAlignments(bf, param = param))) {
    r_total <- .ervmap_3_filters(r, tags_df)
    rm(r)
    overlap <- summarizeOverlaps(features = ervpar@annotations, reads = r_total,
                                 ignore.strand = ervpar@ignoreStrand, mode = Union,
                                 inter.feature = FALSE, singleEnd = TRUE)
    rm(r_total)
    if (is.null(cnt)) {
      cnt <- assay(overlap)
    } else {
      cnt <- cnt + assay(overlap)
    }
  }
  close(bf)
  
  if(tags_df$NH || tags_df$XS) {
    ## If the NH and XS are not provided, both unique reads and multimapping
    ## reads have already been filtered and counted in the previous step. 
    
    sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                           isSupplementaryAlignment = FALSE, isDuplicate = FALSE,
                           isNotPassingQualityControls = FALSE)
    if (!tags_df$NH & tags_df$XS) {
      # If NH is not available but XS is, reads with XS = 0 are considered as unique
      param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"),
                            tagFilter = list(XS = 0))
    } else {
      param <- ScanBamParam(flag = sbflags, tag = c("nM","NM","AS","NH","XS"),
                            tagFilter = list(NH = 1)) #Read only unique reads
    }
    open(bf)
    while (length(r <- readGAlignments(bf, param = param))) {
      if (ervpar@filterUniqReads) {
        r <- .ervmap_2_filters(r, tags_df)
      }
      overlap <- summarizeOverlaps(features = ervpar@annotations, reads = r,
                                   ignore.strand = ervpar@ignoreStrand,
                                   mode = Union, inter.feature = FALSE,
                                   singleEnd = TRUE)
      rm(r)
      cnt <- cnt + assay(overlap)
    }
  }
  return(cnt)
}



#' @importFrom Rsamtools scanBamFlag ScanBamParam yieldSize
#' @importFrom GenomicAlignments readGAlignments summarizeOverlaps last first
#' @importFrom SummarizedExperiment SummarizedExperiment
.qtex_ervmap_pairedend <- function(bf, ervpar) {
  tags_df <- .get_tags_in_BAM_pairedend(bf, ervpar)
  yieldSize(bf) <- 100000
  asMates(bf) <- TRUE
  sbflags <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE,
                         isDuplicate = FALSE, isSupplementaryAlignment = FALSE,
                         isNotPassingQualityControls = FALSE)
  
  if (!tags_df$NH & !tags_df$XS) {
    if (!ervpar@filterUniqReads) {
      stop("Error in 'filterUniqReads = FALSE': Neither the NH tag nor the XS tag (from BWA) are provided in the BAM files. Unique reads cannot be differentiated from multi-mapping reads, therefore unique reads must be also filtered.")
    }
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"))
  } else if (tags_df$NH) {
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"),
                          tagFilter = list(NH = c(2:10000)))
  } else if (tags_df$XS) {
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"),
                          tagFilter = list(XS = c(1:10000)))
  }
  if (!tags_df$NM & !tags_df$nM) {
    stop("Neither the NM tag nor the nM tag are available in the BAM file. At least one of them is needed to apply the 2nd filter.")
  }
  
  cnt <- NULL
  open(bf)
  while (length(r <- readGAlignmentPairs(bf, param = param, strandMode = ervpar@strandMode))) {
    if (ervpar@fragments) {
      r_first_total <- .ervmap_3_filters(first(r), tags_df)
      r_last_total <- .ervmap_3_filters(last(r), tags_df)
      r_total <- c(r_first_total, r_last_total)
      rm(r, r_first_total, r_last_total)
      overlap <- summarizeOverlaps(features = ervpar@annotations, reads = r_total,
                                   ignore.strand = ervpar@ignoreStrand,
                                   mode = Union, inter.feature = FALSE,
                                   fragments = TRUE)
    } else {
      r_total <- .ervmap_3_filters_pairedend(r, tags_df)
      overlap <- summarizeOverlaps(features = ervpar@annotations, reads = r_total,
                                   ignore.strand = ervpar@ignoreStrand,
                                   mode = Union, inter.feature = FALSE,
                                   singleEnd = FALSE, fragments = FALSE)
    }
    rm(r_total)
    if (is.null(cnt)) {
      cnt <- assay(overlap)
    } else {
      cnt <- cnt + assay(overlap)
    }
  }
  close(bf)
  
  if(tags_df$NH || tags_df$XS) {
    ## If the NH and XS are not provided, both unique reads and multimapping
    ## reads have already been filtered and counted in the previous step. 
    sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                           isSupplementaryAlignment = FALSE, 
                           isDuplicate = FALSE, isProperPair = TRUE, 
                           isNotPassingQualityControls = FALSE)
    if (!tags_df$NH & tags_df$XS) {
      # If NH is not available but XS is, reads with XS = 0 are considered as unique
      param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"),
                            tagFilter = list(XS = 0))
    } else {
      param <- ScanBamParam(flag = sbflags, tag = c("nM","NM","AS","NH","XS"),
                            tagFilter = list(NH = 1)) #Read only unique reads
    }
    
    open(bf)
    while (length(r <- readGAlignmentPairs(bf, param = param, strandMode = ervpar@strandMode))) {
      if (ervpar@fragments) {
        if (ervpar@filterUniqReads) {
          r_first_total <- .ervmap_2_filters(first(r), tags_df)
          r_last_total <- .ervmap_2_filters(last(r), tags_df)
          r_total <- c(r_first_total, r_last_total)
          rm(r_first_total, r_last_total)
        } else {
          r_total <- c(first(r), last(r))
        }
        rm(r)
        overlap <- summarizeOverlaps(features = ervpar@annotations, reads = r_total,
                                     ignore.strand = ervpar@ignoreStrand,
                                     mode = Union, inter.feature = FALSE,
                                     fragments = TRUE)
      } else {
        if (ervpar@filterUniqReads) {
          r_total <- .ervmap_2_filters_pairedend(r, tags_df)
        } else {
          r_total <- r
        }
        rm(r)
        overlap <- summarizeOverlaps(features = ervpar@annotations, reads = r_total,
                                     ignore.strand = ervpar@ignoreStrand,
                                     mode = Union, inter.feature = FALSE,
                                     singleEnd = FALSE, fragments = FALSE)
      }
      rm(r_total)
      cnt <- cnt + assay(overlap)
    }
  }
  return(cnt)
}



## Function to know which tags are present in the BAM file (single-end).
## Input is a BamFile object. Returns a data.frame with tags as columns and the 
## BAM file as row, indicating if each tag is present in the BAM file (TRUE) or 
## not (FALSE).

#' @importFrom Rsamtools scanBamFlag ScanBamParam yieldSize
.get_tags_in_BAM_singleend <- function(bf) {
  
  yieldSize(bf) <- 1000
  sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                         isDuplicate = FALSE,
                         isNotPassingQualityControls = FALSE,
                         isSupplementaryAlignment = FALSE)
  
  param <- ScanBamParam(flag = sbflags,
                        tag = c("nM", "NM", "AS", "NH", "XS"))
  
  tags <- param@tag
  tags_df <- data.frame(matrix(nrow = length(bf), ncol = length(tags)))
  colnames(tags_df) <- tags
  
  open(bf)
  r_test <- readGAlignments(bf, param = param)
  close(bf)
  
  # Testing, for each sample, if the files contain the different tags
  tags_df[, "NH"] <- ifelse(all(is.na(mcols(r_test)$NH)), FALSE, TRUE)
  tags_df[, "NM"] <- ifelse(all(is.na(mcols(r_test)$NM)), FALSE, TRUE)
  tags_df[, "nM"] <- ifelse(all(is.na(mcols(r_test)$nM)), FALSE, TRUE)
  tags_df[, "AS"] <- ifelse(all(is.na(mcols(r_test)$AS)), FALSE, TRUE)
  tags_df[, "XS"] <- ifelse(all(grepl("\\d", (mcols(r_test)$XS))), TRUE, FALSE)
  
  return(tags_df)
}



## Function to know which tags are present in the BAM file (paired-end).
## Input is a BamFile object. Returns a data.frame with tags as columns and the 
## BAM file as row, indicating if each tag is present in the BAM file (TRUE) or 
## not (FALSE).

#' @importFrom Rsamtools scanBamFlag ScanBamParam yieldSize asMates
.get_tags_in_BAM_pairedend <- function(bf, ervpar) {
  yieldSize(bf) <- 1000
  asMates(bf) <- TRUE
  sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                         isProperPair = TRUE,
                         isDuplicate = FALSE,
                         isNotPassingQualityControls = FALSE,
                         isSupplementaryAlignment = FALSE)
  
  param <- ScanBamParam(flag = sbflags,
                        tag = c("nM", "NM", "AS", "NH", "XS"))
  
  tags <- param@tag
  tags_df <- data.frame(matrix(nrow = length(bf), ncol = length(tags)))
  colnames(tags_df) <- tags
  
  open(bf)
  r_test <- readGAlignmentPairs(bf, param = param, strandMode = ervpar@strandMode)
  close(bf)
  
  # Testing, for each sample, if the files contain the different tags
  tags_df[, "NH"] <- ifelse(all(is.na(mcols(first(r_test))$NH)), FALSE, TRUE)
  tags_df[, "NM"] <- ifelse(all(is.na(mcols(first(r_test))$NM)), FALSE, TRUE)
  tags_df[, "nM"] <- ifelse(all(is.na(mcols(first(r_test))$nM)), FALSE, TRUE)
  tags_df[, "AS"] <- ifelse(all(is.na(mcols(first(r_test))$AS)), FALSE, TRUE)
  tags_df[, "XS"] <- ifelse(all(grepl("\\d", (mcols(first(r_test))$XS))), TRUE, FALSE)
  
  return(tags_df)
}



## Function to apply the 3 filters from ERVmap to a GAlignments object
## Returns a filtered GAlignments object.

#' @importFrom GenomicAlignments explodeCigarOpLengths qwidth cigar
.ervmap_3_filters <- function(r, tags_df) {
  # The 1st filter is always applied
  cigar_out <- explodeCigarOpLengths(cigar(r), ops = c("H","S"))
  SH_clipping <- lapply(cigar_out, function(cig) sum(unlist(cig)))
  
  # The 2nd filter is always applied (when the NM tag is not available the nM tag is used instead)
  if (tags_df$NM) {
    NM_filter <- (mcols(r)$NM / qwidth(r)) < 0.02
  } else {
    NM_filter <- (mcols(r)$nM / qwidth(r)) < 0.02
  }
  
  # The 3rd filter cannot be applied when neither the XS tag is not available
  if (tags_df$XS) {
    AS_XS_filter <- (mcols(r)$AS - mcols(r)$XS) >= 5
    AS_XS_filter[is.na(AS_XS_filter)] <- FALSE
  } else {
    AS_XS_filter <- rep(TRUE, length(r))
    warning("The XS tag (from BWA) is not provided in the BAM files. This tag is needed to obtain the suboptimal alignment score. The 3rd filter (based on the suboptimal alignment score) will not be applied.")
  }
  
  # Filtering alignments considering all filters
  r_to_keep <- ((unlist(SH_clipping) / qwidth(r)) < 0.02) & NM_filter & AS_XS_filter
  r_total <- r[r_to_keep]
  return(r_total)
}

## Function to apply the first 2 filters from ERVmap to a GAlignments object
## Returns a filtered GAlignments object.

#' @importFrom GenomicAlignments explodeCigarOpLengths qwidth cigar
.ervmap_2_filters <- function(r, tags_df) {
  cigar_out <- explodeCigarOpLengths(cigar(r), ops = c("H","S"))
  SH_clipping <- lapply(cigar_out, function(cig) sum(unlist(cig)))
  if (tags_df$NM) {
    NM_filter <- (mcols(r)$NM / qwidth(r)) < 0.02
  } else if (tags_df$nM) {
    NM_filter <- (mcols(r)$nM / qwidth(r)) < 0.02
  } else {
    stop("Error: Neither the NM tag nor the nM tag are available. At least one of them is needed to apply the 2nd filter.")
  }
  r_to_keep <- ((unlist(SH_clipping) / qwidth(r)) < 0.02) & NM_filter
  r_uniq <- r[r_to_keep]
  return(r_uniq)
}


## Function to apply the 3 filters from ERVmap to a GAlignmentPairs object
## Returns a filtered GAlignmentPairs object.

#' @importFrom GenomicAlignments explodeCigarOpLengths qwidth first last
.ervmap_3_filters_pairedend <- function(r, tags_df) {
  # The first filter is always applied
  cigar_first <- explodeCigarOpLengths(cigar(first(r)), ops = c("H","S"))
  cigar_last <- explodeCigarOpLengths(cigar(last(r)), ops = c("H","S"))
  
  SH_clipping_first <- lapply(cigar_first, function(x) sum(unlist(x)))
  SH_clipping_last <- lapply(cigar_last, function(x) sum(unlist(x)))
  
  # The second filter is always applied but when the NM tag is not
  # available the nM tag can be used instead.
  if (tags_df$NM) {
    NM_filter <- (mcols(first(r))$NM / qwidth(first(r))) < 0.02 & (mcols(last(r))$NM / qwidth(last(r))) < 0.02
    
  } else if (tags_df$nM) {
    NM_filter <- (mcols(first(r))$nM / qwidth(first(r))) < 0.02 & (mcols(last(r))$nM / qwidth(last(r))) < 0.02
  }
  
  # The third filter (AS - XS >= 5) cannot be applied the XS tag is not available
  if (!tags_df$XS) {
    AS_XS_filter <- rep(TRUE, length(r))
  } else {
    AS_XS_filter <- (mcols(first(r))$AS - mcols(first(r))$XS) >= 5 & (mcols(last(r))$AS - mcols(last(r))$XS) >= 5
    AS_XS_filter[is.na(AS_XS_filter)] <- FALSE
  }
  
  # Filtering alignments considering all filters
  r_to_keep <- ((unlist(SH_clipping_first) / qwidth(first(r))) < 0.02 & 
                  (unlist(SH_clipping_last) / qwidth(last(r))) < 0.02 ) & NM_filter & AS_XS_filter
  return(r[r_to_keep])
}



## Function to apply the 2 filters from ERVmap to a GAlignmentPairs object
## Returns a filtered GAlignmentPairs object.

#' @importFrom GenomicAlignments explodeCigarOpLengths qwidth first last
.ervmap_2_filters_pairedend <- function(r, tags_df) {
  # The first filter is always applied
  cigar_first <- explodeCigarOpLengths(cigar(first(r)), ops = c("H","S"))
  cigar_last <- explodeCigarOpLengths(cigar(last(r)), ops = c("H","S"))
  
  SH_clipping_first <- lapply(cigar_first, function(x) sum(unlist(x)))
  SH_clipping_last <- lapply(cigar_last, function(x) sum(unlist(x)))
  
  # The second filter is always applied but when the NM tag is not
  # available the nM tag can be used instead.
  if (tags_df$NM) {
    NM_filter <- (mcols(first(r))$NM / qwidth(first(r))) < 0.02 & (mcols(last(r))$NM / qwidth(last(r))) < 0.02
    
  } else if (tags_df$nM) {
    NM_filter <- (mcols(first(r))$nM / qwidth(first(r))) < 0.02 & (mcols(last(r))$nM / qwidth(last(r))) < 0.02
  }
  
  # Filtering alignments considering all filters
  r_to_keep <- ((unlist(SH_clipping_first) / qwidth(first(r))) < 0.02 & 
                  (unlist(SH_clipping_last) / qwidth(last(r))) < 0.02 ) & NM_filter
  return(r[r_to_keep])
}
