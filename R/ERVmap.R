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
#' empar <- ERVmapParam(bamfiles, annot)
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
              cntmat <- .qtex_ervmap_singleend(x)
            else
              cntmat <- .qtex_ervmap_pairedend(x)

            colnames(cntmat) <- gsub(".bam$", "", colnames(cntmat))
            colData <- DataFrame(row.names=colnames(cntmat))
            if (!is.null(phenodata)) {
              colData <- phenodata
              colnames(cnt) <- rownames(colData)
            }

            SummarizedExperiment(assays=list(counts=cntmat),
                                 rowRanges=x@annotations,
                                 colData=colData)
          })


#' @importFrom Rsamtools scanBamFlag ScanBamParam
.qtex_ervmap_singleend <- function(empar) {
  ## REMOVEME !!
  ## temporary matrix with NAs until the result is properly calculated
  cntmat <- matrix(NA, nrow=length(empar@annotations), ncol=length(empar@bfl),
                   dimnames=list(names(empar@annotations), names(empar@bfl)))

  cntmat
}

#' @importFrom Rsamtools scanBamFlag ScanBamParam
.qtex_ervmap_pairedend <- function(empar) {
  ## REMOVEME !!
  ## temporary matrix with NAs until the result is properly calculated
  cntmat <- matrix(NA, nrow=length(empar@annotations), ncol=length(empar@bfl),
                   dimnames=list(names(empar@annotations), names(empar@bfl)))

  cntmat
}
