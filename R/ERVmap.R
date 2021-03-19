#' Build an ERVmap parameter object
#'
#' Build an object of the class \code{ERVmapParam}
#'
#' @param bfl A \code{BamFile} or \code{BamFileList} object, or a character
#' string vector of BAM filenames. In order to apply all three filters present
#' in the ERVmap algorithm, reads must be aligned using Burrows-Wheeler Aligner.
#' If a different aligner was used, the 3rd filter (AS - XS >= 5) from the 
#' ERVmap algorithm is not applied.
#'
#' @param teFeatures A \code{GRanges} or \code{GRangesList} object with the
#' TE annotated features to be quantified. Elements in this object should have
#' names, which will be used as a grouping factor for genomic ranges forming a
#' common locus, unless other metadata column names are specified in the
#' \code{aggregateby} parameter.
#'
#' @param aggregateby Character vector with column names in the annotation
#' to be used to aggregate quantifications. By default, this is an empty vector,
#' which means that the names of the input \code{GRanges} or \code{GRangesList}
#' object given in the \code{annotations} parameter will be used to aggregate
#' quantifications.
#'
#' @param geneFeatures A \code{GRanges} or \code{GRangesList} object with the
#' gene annotated features to be quantified. Only unique counts are used for
#' quantifying gene features given in this parameter.
#'
#' @param singleEnd (Default TRUE) Logical value indicating if reads are single
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
#' @param fragments (Default not \code{singleEnd}) A logical; applied to
#' paired-end data only. When \code{fragments=TRUE} (default), the read-counting
#' method in the original ERVmap algorithm will be applied, by which each mate
#' of a paired-end read is counted once, and therefore two mates mapping to the
#' same element result in adding up a count value of two. When
#' \code{fragments=FALSE}, if the two mates of a paired-end read map to the same
#' element, they are counted as a single hit and singletons, reads with unmapped
#' pairs and other fragments are not counted.
#'
#' @param filterUniqReads (Default TRUE) Logical value indicating whether to apply
#' the alignment filters to unique reads (TRUE) or not (FALSE). These filters,
#' which are always applied to multi-mapping reads, can be optional for unique
#' reads, only if the NH tag is present in the BAM file. If 
#' \code{filterUniqReads = TRUE} (equivalent to the original approach proposed
#' by ERVmap authors), the unique reads not passing one or more filters 
#' from the ERVmap pipeline will be discarded to compute TEs expression.
#'
#' @details
#' This is the constructor function for objects of the class
#' \code{ERVmapParam-class}. This type of object is the input to the
#' function \code{\link{qtex}()} for quantifying expression of transposable
#' elements using the ERVmap method
#' \href{https://doi.org/10.1073/pnas.1814589115}{Tokuyama et al. (2018)}. The
#' ERVmap algorithm processes reads following conservative filtering criteria
#' to provide reliable raw count data for each TE.
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
ERVmapParam <- function(bfl, teFeatures, aggregateby=character(0),
                        geneFeatures=NA,
                        singleEnd=TRUE,
                        ignoreStrand=TRUE,
                        strandMode=1L,
                        fragments=!singleEnd,
                        filterUniqReads=FALSE) {

  bfl <- .checkBamFileListArgs(bfl, singleEnd, fragments)

  features <- .processFeatures(teFeatures, deparse(substitute(teFeatures)),
                               geneFeatures, deparse(substitute(geneFeatures)),
                               aggregateby)
  
  new("ERVmapParam", bfl=bfl, annotations=features, aggregateby=aggregateby,
      singleEnd=singleEnd, ignoreStrand=ignoreStrand,
      strandMode=as.integer(strandMode), fragments=fragments,
      filterUniqReads=filterUniqReads)
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
#' @export
#' @aliases qtex
#' @aliases qtex,ERVmapParam-method
#' @rdname qtex
setMethod("qtex", "ERVmapParam",
          function(x, phenodata=NULL, mode=ovUnion, usematrix=FALSE,
                   BPPARAM=SerialParam(progressbar=TRUE)) {
            .checkPhenodata(phenodata, length(x@bfl))

            if (usematrix)
              cnt <- bplapply(x@bfl, .qtex_ervmap_matrix, empar=x, BPPARAM=BPPARAM)
            else {
              if (x@singleEnd)
                cnt <- bplapply(x@bfl, .qtex_ervmap_singleend, ervpar=x, BPPARAM=BPPARAM)
              else
                cnt <- bplapply(x@bfl, .qtex_ervmap_pairedend, ervpar=x, BPPARAM=BPPARAM)
            }
            cnt <- do.call("cbind", cnt)
            colData <- .createColumnData(cnt, phenodata)
            colnames(cnt) <- rownames(colData)

            features <- .consolidateFeatures(x, rownames(cnt))

            SummarizedExperiment(assays=list(counts=cnt),
                                 rowRanges=features,
                                 colData=colData)
          })

#' @importFrom Rsamtools scanBamFlag ScanBamParam yieldSize asMates
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom SummarizedExperiment SummarizedExperiment assay
.qtex_ervmap_singleend <- function(bf, ervpar) {
  tags_df <- .get_tags_in_BAM_singleend(bf)
  yieldSize(bf) <- 100000
  sbflags <- scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE,
                         isNotPassingQualityControls = FALSE,
                         isSupplementaryAlignment = FALSE,
                         isSecondaryAlignment = FALSE)
  if (!tags_df$NH) {
    if (!ervpar@filterUniqReads) {
      stop("Error in 'filterUniqReads = FALSE': the NH tag is not provided in the BAM files. Unique reads cannot be differentiated from multi-mapping reads, therefore unique reads must be also filtered.")
    }
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"))
  } else {
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"),
                          tagFilter = list(NH = c(2:10000)))
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
  
  if(tags_df$NH) {
    ## If NH is not provided, both unique reads and multimapping
    ## reads have already been filtered and counted in the previous step. 
    sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                           isSupplementaryAlignment = FALSE, isDuplicate = FALSE,
                           isNotPassingQualityControls = FALSE,
                           isSecondaryAlignment = FALSE)
    param <- ScanBamParam(flag = sbflags, tag = c("nM","NM","AS","NH","XS"),
                          tagFilter = list(NH = 1)) #Read only unique reads
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
#' @importFrom SummarizedExperiment SummarizedExperiment assay
.qtex_ervmap_pairedend <- function(bf, ervpar) {
  tags_df <- .get_tags_in_BAM_pairedend(bf, ervpar)
  yieldSize(bf) <- 100000
  asMates(bf) <- TRUE
  sbflags <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE,
                         isDuplicate = FALSE, isSupplementaryAlignment = FALSE,
                         isNotPassingQualityControls = FALSE,
                         isSecondaryAlignment = FALSE)
  
  if (!tags_df$NH) {
    if (!ervpar@filterUniqReads) {
      stop("Error in 'filterUniqReads = FALSE': the NH tag is not provided in the BAM files. Unique reads cannot be differentiated from multi-mapping reads, therefore unique reads must be also filtered.")
    }
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"))
  } else {
    param <- ScanBamParam(flag = sbflags, tag = c("nM", "NM", "AS", "NH", "XS"),
                          tagFilter = list(NH = c(2:10000)))
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
      overlap <- summarizeOverlaps(features = ervpar@annotations,
                                   reads = r_total, 
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
  
  if(tags_df$NH) {
    ## If the NH is not provided, both unique reads and multimapping
    ## reads have already been filtered and counted in the previous step. 
    sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                           isSupplementaryAlignment = FALSE, 
                           isDuplicate = FALSE, isProperPair = TRUE, 
                           isNotPassingQualityControls = FALSE,
                           isSecondaryAlignment = FALSE)
    param <- ScanBamParam(flag = sbflags, tag = c("nM","NM","AS","NH","XS"),
                          tagFilter = list(NH = 1)) #Read only unique reads
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
#' @importFrom S4Vectors mcols
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
#' @importFrom S4Vectors mcols
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
#' @importFrom S4Vectors mcols
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
#' @importFrom S4Vectors mcols
.ervmap_2_filters <- function(r, tags_df) {
  cigar_out <- explodeCigarOpLengths(cigar(r), ops = c("H","S"))
  SH_clipping <- lapply(cigar_out, function(cig) sum(unlist(cig)))
  if (tags_df$NM) {
    NM_filter <- (mcols(r)$NM / qwidth(r)) < 0.02
  } else {
    NM_filter <- (mcols(r)$nM / qwidth(r)) < 0.02
  }
  r_to_keep <- ((unlist(SH_clipping) / qwidth(r)) < 0.02) & NM_filter
  r_uniq <- r[r_to_keep]
  return(r_uniq)
}


## Function to apply the 3 filters from ERVmap to a GAlignmentPairs object
## Returns a filtered GAlignmentPairs object.

#' @importFrom GenomicAlignments explodeCigarOpLengths qwidth first last
#' @importFrom S4Vectors mcols
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
  } else {
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
#' @importFrom S4Vectors mcols
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
  } else {
    NM_filter <- (mcols(first(r))$nM / qwidth(first(r))) < 0.02 & (mcols(last(r))$nM / qwidth(last(r))) < 0.02
  }
  
  # Filtering alignments considering all filters
  r_to_keep <- ((unlist(SH_clipping_first) / qwidth(first(r))) < 0.02 & 
                  (unlist(SH_clipping_last) / qwidth(last(r))) < 0.02 ) & NM_filter
  return(r[r_to_keep])
}

#' @importFrom Rsamtools scanBamFlag ScanBamParam yieldSize
#' @importFrom Matrix rowSums
.qtex_ervmap_matrix <- function(bf, empar, mode) {
  
  mode=match.fun(mode)

  readfun <- .getReadFunction(empar@singleEnd, empar@fragments)

  sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                         isDuplicate=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isSupplementaryAlignment=FALSE)
  param <- ScanBamParam(flag=sbflags, what="flag", tag=c("AS", "NM"))

  ov <- Hits(nLnode=0, nRnode=length(empar@annotations), sort.by.query=TRUE)
  alnreadids <- character(0)
  salnmask <- logical(0)
  alnAS <- alnNM <- integer(0)

  yieldSize(bf) <- 100000
  open(bf)
  while (length(alnreads <- readfun(bf, param=param, use.names=TRUE))) {
    aqw <- .getAlignmentQueryWidth(alnreads)
    anm <- .getAlignmentMismatches(alnreads)
    asc <- .getAlignmentSumClipping(alnreads)
    ## Tokuyama et al. (2018) pg. 12571. reads/fragments are kept
    ## if: (1) the ratio of sum of hard and soft clipping to the
    ## sequence read length is < 0.02; (2) the ratio of the edit
    ## distance to the sequence read length is < 0.02
    mask <- ((asc / aqw) < 0.02) & ((anm / aqw) < 0.02)
    alnreads <- alnreads[mask, ]
    anm <- anm[mask]

    alnreadids <- c(alnreadids, names(alnreads))
    salnmask <- c(salnmask, .secondaryAlignmentMask(alnreads))
    alnAS <- c(alnAS, .getAlignmentScore(alnreads))
    alnNM <- c(alnNM, .getAlignmentMismatches(alnreads))

    thisov <- mode(alnreads, empar@annotations, ignoreStrand=empar@ignoreStrand)
    ov <- .appendHits(ov, thisov)
  }
  close(bf)

  ## fetch all different read identifiers from the overlapping alignments
  readids <- unique(alnreadids[queryHits(ov)])

  ## fetch all different transcripts from the overlapping alignments
  tx_idx <- sort(unique(subjectHits(ov)))

  ## build a matrix representation of the overlapping alignments
  ovalnmat <- .buildOvAlignmentsMatrix(ov, alnreadids, readids, tx_idx)

  ## Tokuyama et al. (2018) pg. 12571. reads/fragments are kept
  ## if the difference between the alignment score from BWA (AS)
  ## and the suboptimal alignment score from BWA (field XS) >= 5,
  ## equivalent to second best match has one or more mismatches
  ## than the best match. here we pick the best suboptimal alignment
  ## as the secondary alignment with highest AS and compare the
  ## number of mismatches (NM) to the primary alignment for that read.
  ## if NM-primary-alignment - NM-best-suboptimal-alignment >= 1,
  ## then the read/fragment is kept.
  if (any(salnmask)) {
    salnmat <- .buildSecondaryAlignmentsMatrix(ov, alnreadids, salnmask, readids, tx_idx)
    asmat <- .buildAlignmentScoresMatrix(ov, alnreadids, alnAS, readids, tx_idx)
    nmmat <- .buildAlignmentScoresMatrix(ov, alnreadids, alnNM, readids, tx_idx)
    ## remove matrices rownames to save memory? (they're 'readids')

    nmprimaryaln <- rowMeans(nmmat * salnmat)
    whsalnmaxas <- max.col(salnmat * asmat)
    nmbestsecondaryaln <- nmmat[, whsalnmaxas]

    mask <- nmprimaryaln - nmbestsecondaryaln >= 1
    if (!any(mask))
      stop("no read passes filters.")

    ovalnmat <- ovalnmat[mask, ]
  }

  rm(alnreadids)
  gc()

  cntvec <- rep(0, length(empar@annotations))
  cntvec[tx_idx] <- colSums(ovalnmat)
  names(cntvec) <- names(empar@annotations)

  ## aggregate quantifications if necessary
  if (length(empar@aggregateby) > 0) {
    f <- .factoraggregateby(empar@annotations, empar@aggregateby)
    stopifnot(length(f) == length(cntvec)) ## QC
    cntvec <- tapply(cntvec, f, sum, na.rm=TRUE)
  }

  cntvec
}

#' @importFrom S4Vectors mcols first second
#' @importFrom Rsamtools bamFlagTest
.secondaryAlignmentMask <- function(aln) {
  mask <- NULL
  if (is(aln, "GAlignments"))
    mask <- bamFlagTest(mcols(aln)$flag, "isSecondaryAlignment")
  else if (is(aln, "GAlignmentPairs")) {
    mask <- bamFlagTest(mcols(first(aln))$flag, "isSecondaryAlignment") |
            bamFlagTest(mcols(second(aln))$flag, "isSecondaryAlignment")
  } else if (is(aln, "GAlignmentsList")) {
    maskl <- relist(bamFlagTest(mcols(unlist(aln, use.names=FALSE))$flag,
                                "isSecondaryAlignment"), aln)
    names(maskl) <- NULL
    ## any secondary alignment makes all the grouped alignments secondary
    mask <- sum(maskl) > 0
  } else
    stop(sprintf(".secondaryAlignmentMask: wrong class %s\n", class(aln)))

  mask
}

#' @importFrom S4Vectors mcols first second
.getAlignmentScore <- function(aln) {
  if (is(aln, "GAlignments"))
    ascore <- as.integer(mcols(aln)$AS)
  else if (is(aln, "GAlignmentPairs")) {
    ## take the minimum score from both mates
    ascore <- pmin.int(as.integer(mcols(first(aln))$AS), mcols(second(aln))$AS)
  } else if (is(aln, "GAlignmentsList")) {
    ascorel <- relist(as.integer(mcols(unlist(aln, use.names=FALSE))$AS), aln)
    names(ascorel) <- NULL
    ## any secondary alignment makes all the grouped alignments secondary
    ascore <- min(ascorel)
  } else
    stop(sprintf(".getAlignmentScore: wrong class %s\n", class(aln)))

  as.integer(ascore)
}

.fetchNMtag <- function(aln) {
  nmtag <- "NM"
  if (is(aln, "GAlignments")) {
    if (is.null(mcols(aln)[[nmtag]]) || all(is.na(mcols(aln)[[nmtag]]))) {
      nmtag <- "nM"
      if (is.null(mcols(aln)[[nmtag]]) || all(is.na(mcols(aln)[[nmtag]])))
        stop("no NM or nM tag in BAM file.")
    }
  } else if (is(aln, "GAlignmentPairs")) {
    if (is.null(mcols(first(aln))[[nmtag]]) || all(is.na(mcols(first(aln))[[nmtag]]))) {
      nmtag <- "nM"
      if (is.null(mcols(first(aln))[[nmtag]]) || all(is.na(mcols(first(aln))[[nmtag]])))
        stop("no NM or nM tag in BAM file.")
    }
  } else if (is(aln, "GAlignmentsList")) {
    if (is.null(mcols(unlist(aln))[[nmtag]]) || all(is.na(mcols(unlist(aln))[[nmtag]]))) {
      nmtag <- "nM"
      if (is.null(mcols(unlist(aln))[[nmtag]]) || all(is.na(mcols(unlist(aln))[[nmtag]])))
        stop("no NM or nM tag in BAM file.")
    }
  } else
    stop(sprintf(".getAlignmentMismatches: wrong class %s\n", class(aln)))

  nmtag
}

#' @importFrom S4Vectors mcols first second
.getAlignmentMismatches <- function(aln) {
  mm <- NULL
  nmtag <- .fetchNMtag(aln)

  if (is(aln, "GAlignments"))
    mm <- as.integer(mcols(aln)[[nmtag]])
  else if (is(aln, "GAlignmentPairs")) {
    ## take the celing of the mean of mismatches from both mates
    mm1 <- as.integer(mcols(first(aln))[[nmtag]])
    mm2 <- as.integer(mcols(second(aln))[[nmtag]])
    mm1[is.na(mm1)] <- mm2[is.na(mm1)]
    mm2[is.na(mm2)] <- mm1[is.na(mm2)]
    mm <- ceiling((mm1 + mm2) / 2)
  } else if (is(aln, "GAlignmentsList")) {
    mml <- relist(as.integer(mcols(unlist(aln, use.names=FALSE))[[nmtag]]), aln)
    names(mml) <- NULL
    len <- length(mml)
    ## take the celing of the mean of mismatches from the alignment group
    mm <- ceiling(sum(mml) / len)
  } else
    stop(sprintf(".getAlignmentMismatches: wrong class %s\n", class(aln)))

  as.integer(mm)
}

#' @importFrom S4Vectors mcols first second
#' @importFrom GenomicAlignments qwidth
.getAlignmentQueryWidth <- function(aln) {
  qw <- NULL
  if (is(aln, "GAlignments"))
    qw <- qwidth(aln)
  else if (is(aln, "GAlignmentPairs")) {
    ## take the ceiling of the mean of query widths from both mates
    qw1 <- qwidth(first(aln))
    qw2 <- qwidth(second(aln))
    qw1[is.na(qw1)] <- qw2[is.na(qw1)]
    qw2[is.na(qw2)] <- qw1[is.na(qw2)]
    qw <- ceiling((qw1 + qw2) / 2)
  } else if (is(aln, "GAlignmentsList")) {
    qwl <- relist(qwidth(unlist(aln, use.names=FALSE)), aln)
    names(qwl) <- NULL
    len <- lengths(qwl)
    ## take the ceiling of the mean of query widths from the alignment group
    qw <- ceiling(sum(qwl) / len)
  } else
    stop(sprintf(".getAlignmentQueryWidth: wrong class %s\n", class(aln)))

  as.integer(qw)
}

#' @importFrom S4Vectors mcols first second List sum
#' @importFrom GenomicAlignments cigar explodeCigarOpLengths
.getAlignmentSumClipping <- function(aln) {
  sc <- NULL
  if (is(aln, "GAlignments"))
    sc <- unlist(sum(List(explodeCigarOpLengths(cigar(aln), ops=c("H", "S")))),
                 use.names=FALSE)
  else if (is(aln, "GAlignmentPairs")) {
    ## take the ceiling of the mean of query widths from both mates
    sc1 <- unlist(sum(List(explodeCigarOpLengths(cigar(first(aln)),
                                                 ops=c("H", "S")))),
                  use.names=FALSE)
    sc2 <- unlist(sum(List(explodeCigarOpLengths(cigar(second(aln)),
                                                 ops=c("H", "S")))),
                  use.names=FALSE)
    sc1[is.na(qw1)] <- sc2[is.na(sc1)]
    sc2[is.na(qw2)] <- sc1[is.na(sc2)]
    sc <- ceiling((sc1 + sc2) / 2)
  } else if (is(aln, "GAlignmentsList")) {
    scl <- List(explodeCigarOpLengths(cigar(unlist(aln, use.names=FALSE)),
                                      ops=c("H", "S")))
    scl <- relist(sum(scl), aln)
    names(scl) <- NULL
    len <- lengths(scl)
    ## take the ceiling of the mean of soft clippings from the alignment group
    sc <- ceiling(sum(scl) / len)
  } else
    stop(sprintf(".getAlignmentSumClipping: wrong class %s\n", class(aln)))

  as.integer(sc)
}

#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix Matrix
.buildAlignmentMismatchesMatrix <- function(ov, arids, alnnm, rids, fidx) {
  nmmat <- Matrix(0L, nrow=length(rids), ncol=length(fidx),
                  dimnames=list(rids, NULL))
  mt1 <- match(arids[queryHits(ov)], rids)
  mt2 <- match(subjectHits(ov), fidx)
  nmmat[cbind(mt1, mt2)] <- alnnm[queryHits(ov)]

  nmmat
}

#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix Matrix
.buildAlignmentScoresMatrix <- function(ov, arids, alnas, rids, fidx) {
  asmat <- Matrix(0L, nrow=length(rids), ncol=length(fidx),
                  dimnames=list(rids, NULL))
  mt1 <- match(arids[queryHits(ov)], rids)
  mt2 <- match(subjectHits(ov), fidx)
  asmat[cbind(mt1, mt2)] <- alnas[queryHits(ov)]

  asmat
}

#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix Matrix
.buildSecondaryAlignmentsMatrix <- function(ov, arids, samask, rids, fidx) {
  samat <- Matrix(FALSE, nrow=length(rids), ncol=length(fidx),
                  dimnames=list(rids, NULL))
  mt1 <- match(arids[queryHits(ov)], rids)
  mt2 <- match(subjectHits(ov), fidx)
  samat[cbind(mt1, mt2)] <- samask[queryHits(ov)]

  samat
}
