#' Build an ERVmap parameter object
#'
#' This function serves the purpose of building an object of the
#' class \code{ERVmapParam-class}.
#'
#' @param bfl A 'BamFile' or 'BamFileList' object, or a character string vector
#' of BAM filenames.
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
#' @return A \code{ERVmapParam-class} object.
#'
#' @examples
#' bamfiles <- list.files(system.file("extdata", package="atena"),
#'                        pattern="*.bam", full.names=TRUE)
#' annot <- ERVmap_ann()
#' empar <- ERVmapParam(bamfiles, annot)
#' empar
#'
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

  new("ERVmapParam", bfl=bfl, annotations=annotations, singleEnd=singleEnd,
      ignoreStrand=ignoreStrand, strandMode=as.integer(strandMode),
      filterUniqReads=filterUniqReads, fragments=fragments)
}

#' @importFrom BiocGenerics path
#' @export
#' @aliases path,AtenaParam-method
#' @rdname ERVmapParam-class
setMethod("path", "AtenaParam",
          function(object) {
            path(object@bfl)
          })

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

#' Quantify transposable element expression
#'
#' Giving an \code{AtenaParam} object as input, the
#' \code{qtex()} method quantifies the expression of transposable
#' elements (TEs). The particular algorithm to perform the
#' quantification will be selected depending on the specific
#' sub-class of input \code{AtenaParam} object, see argument
#' \code{x} below.
#'
#' @param x An \code{AtenaParam} object of one of the following
#' subclasses:
#' \itemize{
#'   \item A \code{ERVmapParam} object built using the constructor
#'         function \code{\link{ERVmapParam}()}. This object will
#'         trigger \code{qtex()} to use the algorithm by
#'         Tokuyama et al. (2018).
#'   \item A \code{TelescopeParam} object built using the constructor
#'         function \code{\link{TelescopeParam}()}. This object will
#'         trigger \code{qtex()} to use the algorithm by
#'         Bendall et al. (2019).
#' }
#'
#' @return A \code{SummarizedExperiment} object.
#'
#' @seealso
#' \code{\link{ERVmapParam}}

#' @references
#' Tokuyama M et al. ERVmap analysis reveals genome-wide transcription of human
#' endogenous retroviruses. PNAS, 115(50):12565-12572, 2018.
#' \url{https://doi.org/10.1073/pnas.1814589115}
#'
#' @references
#' Bendall ML et al. Telescope: characterization of the retrotranscriptome by
#' accurate estimation of transposable element expression.
#' PLOS Computational Biology, 15:e1006453, 2019.
#' \url{https://doi.org/10.1371/journal.pcbi.1006453}
#'
#' @export
#' @aliases qtex
#' @aliases qtex,ERVmapParam-method
#' @rdname qtex
setMethod("qtex", "ERVmapParam",
          function(x) {
            cat("ERVmap!!\n")
          })
