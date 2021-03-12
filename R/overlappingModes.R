#' Pre-defined overlapping mode functions
#'
#' The following functions control the way in which overlaps between
#' aligned reads and annotated features are resolved when an aligned
#' read overlaps more than one feature on the same locus:
#'
#' \itemize{
#'   \item \code{ovUnion()}: (default)
#'   \item \code{ovIntersectionStrict()}:
#'   \item \code{ovIntersectionNotEmpty()}:
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
#' @details
#'
#' These functions are given to the \code{mode} parameter of the
#' \code{\link{qtex}()} function and are similar to the functions
#' \code{\link[GenomicAlignments]{Union}()},
#' \code{\link[GenomicAlignments]{IntersectionStrict}()} and
#' \code{\link[GenomicAlignments]{IntersectionNonEmpty}()} from the
#' \code{GenomicAlignments} package, with the difference that instead
#' of returning counts of reads overlapping annotated features, they
#' return the actual overlaps, because the counting is deferred to other
#' algorithms that follow some specific strategy when a read maps to
#' more than one feature. For this same reason, these functions lack
#' the \code{inter.feature} argument found in the corresponding functions
#' from the \code{GenomicAlignments} package.
#'
#' @importFrom GenomicAlignments findOverlaps
#' @importFrom S4Vectors countQueryHits
#' @export
#' @rdname ovFunctions
ovUnion <- function(reads, features, ignoreStrand) {
  ov <- findOverlaps(reads, features, ignore.strand=ignoreStrand)
  ov
}

#' @importFrom GenomicAlignments findOverlaps
#' @importFrom S4Vectors countQueryHits
#' @export
#' @rdname ovFunctions
ovIntersectionStrict <- function(reads, features, ignoreStrand) {
  ov <- findOverlaps(reads, features, type="within",
                     ignore.strand=ignoreStrand)
  ov
}

#' @importFrom GenomicAlignments findOverlaps
#' @importFrom S4Vectors countQueryHits
#' @importFrom GenomicRanges disjoin
.removeSharedRegions <- function(features, ignore.strand=FALSE) {
    if (is(features, "GRanges")) {
        regions <- disjoin(features, ignore.strand=ignore.strand)
    } else if (is(features, "GRangesList")) {
        regions <- disjoin(features@unlistData, ignore.strand=ignore.strand)
    } else {
        stop("internal error")  # should never happen
    }
    ov <- findOverlaps(features, regions,
                       ignore.strand=ignore.strand)
    regions_to_keep <- which(countSubjectHits(ov) == 1L)
    ov <- ov[subjectHits(ov) %in% regions_to_keep]
    unlisted_ans <- regions[subjectHits(ov)]
    ans_partitioning <- as(ov, "PartitioningByEnd")
    relist(unlisted_ans, ans_partitioning)
}

#' @importFrom GenomicAlignments findOverlaps
#' @importFrom S4Vectors countQueryHits
#' @export
#' @rdname ovFunctions
ovIntersectionNonEmpty <- function(reads, features, ignoreStrand) {
  features <- .removedSharedRegions(features, ignore.strand=ignoreStrand)
  ovUnion(reads, features, ignoreStrand=ignoreStrand)
}
