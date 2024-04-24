#' Pre-defined overlapping mode functions
#'
#' The following functions control the way in which overlaps between
#' aligned reads and annotated features are resolved when an aligned
#' read overlaps more than one feature on the same locus:
#'
#' \itemize{
#'   \item \code{ovUnion()}: (default)
#'   \item \code{ovIntersectionStrict()}:
#'   \item User supplied: a function taking the same parameters as the
#'         previous three functions and returning a
#'         \code{\link[S4Vectors]{Hits}} object.
#' }
#'
#' They take the following parameters:
#'
#' @param reads A \code{GAlignments}, \code{GAlignmentList} or a
#'        \code{GAlignmentPairs} object.
#'
#' @param features A \code{GRanges} object with annotated features.
#'
#' @param ignoreStrand (Default FALSE) A logical which defines if the strand
#' should be taken into consideration when computing the overlap between reads
#' and annotated features. When \code{ignoreStrand = FALSE}, an aligned read
#' will be considered to be overlapping an annotated feature as long as they
#' have a non-empty intersecting genomic ranges on the same strand, while when
#' \code{ignoreStrand = TRUE} the strand will not be considered.
#' 
#' @param inter.feature When TRUE, ambiguous alignments (alignments
#' overlapping > 1 features) are removed and not counted. When 
#' \code{inter.feature} is set to FALSE, these ambiguous overlaps are taken
#' into account and addressed differently depending on the TE quantification.
#' 
#'
#' @details These functions are given to the \code{mode} parameter of the
#' \code{\link{qtex}()} function and are similar to the functions
#' \code{\link[GenomicAlignments]{Union}()} and
#' \code{\link[GenomicAlignments]{IntersectionStrict}()} from the
#' \code{GenomicAlignments} package, with the difference that instead
#' of returning counts of reads overlapping annotated features, they
#' return the actual overlaps, because the counting is deferred to other
#' algorithms that follow some specific strategy when a read maps to
#' more than one feature. For this same reason, these functions lack
#' the \code{inter.feature} argument found in the corresponding functions
#' from the \code{GenomicAlignments} package.
#'
#' @return A \code{Hits} object; see the \code{\link[S4Vectors]{Hits-class}}
#' manual page.
#'
#' @examples
#' bamfiles <- list.files(system.file("extdata", package="atena"),
#'                        pattern="*.bam", full.names=TRUE)
#' \dontrun{
#' ## use the following two instructions to fetch annotations, they are here
#' ## commented out to enable running this example quickly when building and
#' ## checking the package
#' rmskat <- annotaTEs(genome="dm6", parsefun=rmskatenaparser,
#'                     strict=FALSE, insert=500)
#' rmskLTR <- getLTRs(rmskat, relLength=0.8, fullLength=TRUE, partial=TRUE,
#'                    otherLTR=TRUE)
#' }
#'
#' ## DO NOT TYPE THIS INSTRUCTION, WHICH JUST LOADS A PRE-COMPUTED ANNOTATION
#' ## YOU SHOULD USE THE INSTRUCTIONS ABOVE TO FETCH ANNOTATIONS
#' rmskLTR <- readRDS(system.file("extdata", "rmskatLTRrlen80flenpartoth.rds",
#'                                package="atena"))
#'
#' ## build a parameter object for Telescope
#' tspar <- TelescopeParam(bfl=bamfiles,
#'                         teFeatures=rmskLTR,
#'                         singleEnd=TRUE,
#'                         ignoreStrand=TRUE)
#'
#' ## quantify expression using the 'ovIntersectionStrict()' mode function
#' tsquant <- qtex(tspar, mode=ovIntersectionStrict)
#'
#' @importFrom GenomicAlignments findOverlaps
#' @importFrom S4Vectors countQueryHits
#' @export
#' @rdname ovFunctions
ovUnion <- function(reads, features, ignoreStrand, inter.feature=TRUE) {
    ov <- findOverlaps(reads, features, ignore.strand=ignoreStrand)
    if (inter.feature) {
      ## Remove ambigous reads.
      reads_to_keep <- which(countQueryHits(ov) == 1L)
      ov <- ov[queryHits(ov) %in% reads_to_keep]
    }
    
    ov
}


#' @importFrom GenomicAlignments findOverlaps
#' @importFrom S4Vectors countQueryHits
#' @export
#' @rdname ovFunctions
ovIntersectionStrict <- function(reads, features, ignoreStrand, 
                                 inter.feature=TRUE) {
    ov <- findOverlaps(reads, features, type="within",
                       ignore.strand=ignoreStrand)
    if (inter.feature) {
      ## Remove ambigous reads.
      reads_to_keep <- which(countQueryHits(ov) == 1L)
      ov <- ov[queryHits(ov) %in% reads_to_keep]
    }
    
    ov
}
