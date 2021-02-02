#' Atena parameter class
#'
#' This is a virtual class from which other classes are derived
#' for storing parameters provided to quantification methods of
#' transposable elements from RNA-seq data.
#'
#' @slot bfl A \linkS4class{BamFileList} object.
#'
#' @slot annotations A \linkS4class{GRanges} object.
#'
#' @importClassesFrom Rsamtools BamFileList
#' @importClassesFrom GenomicRanges GenomicRanges_OR_GenomicRangesList
#'
#' @seealso
#' \code{\link{ERVmapParam-class}}
#' \code{\link{TelescopeParam-class}}
#'
#' @name AtenaParam-class
#' @rdname AtenaParam-class
#' @exportClass AtenaParam
setClass("AtenaParam",
         representation(bfl="BamFileList",
                        annotations="GenomicRanges_OR_GenomicRangesList"))

#' @param object A \linkS4class{AtenaParam} object.
#'
#' @importFrom BiocGenerics path
#' @export
#' @aliases path,AtenaParam-method
#' @rdname AtenaParam-class
setMethod("path", "AtenaParam",
          function(object) {
            path(object@bfl)
          })

#' ERVmap parameter class
#'
#' This is a class for storing parameters provided to the ERVmap algorithm.
#' It is a subclass of the 'AtenaParam-class'.
#'
#' @slot singleEnd (Default FALSE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @slot strandMode (Default 1) Numeric vector which can take values 0, 1 or 2.
#'   The strand mode is a per-object switch on
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   objects that controls the behavior of the strand getter. See
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   class for further detail. If \code{singleEnd = TRUE}, then use either
#'   \code{strandMode = NULL} or do not specify the \code{strandMode} parameter.
#'
#' @slot ignoreStrand (Default TRUE) A logical which defines if the strand
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
#' @slot filterUniqReads (Default TRUE) Logical value indicating whether to apply
#' the alignment filters to unique reads (TRUE) or not (FALSE). These filters,
#' which are always applied to multi-mapping reads, are optional for unique
#' reads. If TRUE, the unique reads not passing one or more filters from the
#' ERVmap pipeline, except for the "AS - XS >= 5" filter, will be discarded to
#' compute TEs expression.
#'
#' @slot fragments (Default TRUE) A logical; applied to paired-end data only.
#' When \code{fragments=TRUE} (default), the read-counting method in the
#' original ERVmap algorithm will be applied, by which each mate of a paired-end
#' read is counted once, and therefore two mates mapping to the same element
#' result in adding up a count value of two. When \code{fragments=FALSE}, if the
#' two mates of a paired-end read map to the same element, they are counted as
#' a single hit and singletons, reads with unmapped pairs and other fragments
#' are not counted.
#' 
#' @references
#' Tokuyama M et al. ERVmap analysis reveals genome-wide transcription of human
#' endogenous retroviruses. PNAS. 2018;115(50):12565-12572. DOI:
#' \url{https://doi.org/10.1073/pnas.1814589115}
#'
#' @name ERVmapParam-class
#' @rdname ERVmapParam-class
#' @exportClass ERVmapParam
setClass("ERVmapParam", contains="AtenaParam",
         representation(singleEnd="logical",
                        ignoreStrand="logical",
                        strandMode="integer",
                        filterUniqReads="logical",
                        fragments="logical"))

#' Telescope parameter class
#'
#' This is a class for storing parameters
#' provided to the Telescope algorithm.
#'
#' @slot basiliskEnv A 'BasiliskEnvironment' object; see the
#' \code{BasiliskEnvironment-class} in the \code{basilisk} package.
#'
#' @slot telescopeVersion 'character' string storing the Telescope version.
#'
#' @slot telescopeOptions 'list' object storing the options for the
#' Telescope algorithm.
#'
#' @references
#' Bendall et al. Telescope: characterization of the retrotranscriptome by
#' accurate estimation of transposable element expression.
#' PLOS Comp. Biol. 2019;15(9):e1006453. DOI:
#' \url{https://doi.org/10.1371/journal.pcbi.1006453}
#'
#' @name TelescopeParam-class
#' @rdname TelescopeParam-class
#' @importClassesFrom basilisk BasiliskEnvironment
#' @exportClass TelescopeParam
setClass("TelescopeParam", contains="AtenaParam",
         representation(basiliskEnv="BasiliskEnvironment",
                        telescopeVersion="character",
                        telescopeOptions="list"))
