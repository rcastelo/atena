#' Obtain the suboptimal alignment score tag
#'
#' \code{getSuboptimalAlignScore} computes the suboptimal alignment score for
#' primary alignments by leveraging the information provided by
#' secondary alignments of multi-mapping reads in binary 'BAM' files. In
#' summary, the secondary alignments of each primary alignment are identified
#' and the alignment score of the best secondary alignment is reported as the
#' suboptimal alignment score. This function is called internally by
#' \code{\link{count_TE}} when the suboptimal alignment score tag is not
#' reported in the provided BAM files, since this information is required by
#' \code{\link{count_TE}} to correctly quantify the expression of transposable
#' elements.
#'
#' This function aims at replacing the 'XS' tag reported by the
#' \href{http://bio-bwa.sourceforge.net/bwa.shtml}{Burrows-Wheeler Aligner}
#' (BWA), so that RNA-sequencing datasets aligned to the reference
#' genome using other software tools different to BWA can also take advantage
#' of this pipeline.
#'
#'
#' @param r \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   object lacking the suboptimal alignment score tag and containing only
#'   alignments corresponding to multi-mapping reads.
#'
#' @param strandMode (Default 1) Numeric vector which can take values 0, 1 or 2. The strand
#'   mode is a per-object switch on
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   objects that controls the behavior of the strand getter. See
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   class for further detail.
#'
#' @return A \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   object with primary alignments of multi-mapping reads displaying the
#'   suboptimal alignment score tag ('XS') for both the first and last mate of
#'   the paired-end reads.
#'
#' @seealso
#' \code{\link[Rsamtools]{BamFile}},
#' \code{\link[Rsamtools]{scanBamFlag}},
#' \code{\link[Rsamtools]{ScanBamParam}},
#' \code{\link[GenomicAlignments]{readGAlignmentPairs}},
#' \code{\link[GenomicAlignments]{GAlignmentPairs-class}},
#' \code{\link{count_TE}}
#'
#' @examples
#' \dontrun{path_bam_files_sorted <- "Documents/RNAseq/BAMfiles_sorted"
#' bam_files <- list.files(path_bam_files_sorted,
#'                         pattern = "\\.bam$",
#'                         full.names = TRUE)
#' bfl <- BamFileList(bam_files, asMates = TRUE)
#'
#' sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
#'                        isProperPair = TRUE,
#'                        isDuplicate = FALSE,
#'                        isNotPassingQualityControls = FALSE)
#'
#' param <- ScanBamParam(flag = sbflags,
#'                       tag = c("NM", "AS", "NH", "XS"),
#'                       what = scanBamWhat(),
#'                       tagFilter = list(NH = c(2:10000)))
#'
#' r <- readGAlignmentPairs(file, use.names=TRUE,
#'                          param=param, strandMode=2)
#'
#' getSuboptimalAlignScore(r, strandMode=2)}
#'
#' @export
#' @import GenomicAlignments
#' @import Rsamtools


getSuboptimalAlignScore <- function(r, strandMode = 1) {

  # -- First, find the suboptimal alignment score for the "first" mates

  if (all(is.na(mcols(first(r))$NH))) {
    stop("The GAlignmentPairs object must have the 'NH' tag")
  }

  r_first <- first(r)

  if (!is.integer(mcols(r_first)$flag))
    stop("The 'flag' column must be an integer vector")

  secondary <- bamFlagAsBitMatrix(mcols(r_first)$flag)[, "isSecondaryAlignment"]

  mcols_r <- mcols(r_first)

  for (n in n:length(r_first)) {

    # if it is a primary alignment (secondary == 0) then:
    if (secondary[n] == 0) {

      # Determining the initial and final positions in which the secondary
      # alignments of a primary alignment will be searched
      initial_pos <- ifelse(n - 100 < 1, 1, n - 100)
      final_pos <- ifelse(n + 100 > length(r_first),  length(r_first), n + 100)
      pos <- initial_pos:final_pos
      pos <- pos[!pos == n]

      # Obtaining the suboptimal alignment score
      secondary <- which(mcols_r[pos, "qname"] == mcols_r[n, "qname"] &
                           as.vector(mcols_r[pos, "seq"]) ==
                           as.vector(mcols_r[n, "seq"]))

      mcols(r_first[n])$XS <- max(mcols_r[secondary, "AS"])

    }

  }

  # -- Second, find the suboptimal alignment score for the "last" mates

  # Selecting alignments corresponding to multi-mapping reads
  if (all(is.na(mcols(last(r))$NH))) {
    stop("The GAlignmentPairs object must have the 'NH' tag")
  }

  r_last <- last(r)

  if (!is.integer(mcols(r_last)$flag))
    stop("The 'flag' column must be an integer vector")

  secondary <- bamFlagAsBitMatrix(mcols(r_last)$flag)[, "isSecondaryAlignment"]

  mcols_r <- mcols(r_last)

  for (n in n:length(r_last)) {

    # if it is a primary alignment (secondary == 0)
    if (secondary[n] == 0) {

      # Determining the initial and final positions in which the secondary
      # alignments of a multi-mapping read will be searched
      initial_pos <- ifelse(n - 100 < 1,  1, n - 100)
      final_pos <- ifelse(n + 100 > length(r_last),  length(r_last), n + 100)
      pos <- initial_pos:final_pos
      pos <- pos[!pos == n]

      # Obtaining the suboptimal alignment score
      secondary <- which(mcols_r[pos, "qname"] == mcols_r[n, "qname"] &
                           as.vector(mcols_r[pos, "seq"]) ==
                           as.vector(mcols_r[n, "seq"]))

      mcols(r_last[n])$XS <- max(mcols_r[secondary, "AS"])

    }

  }

  # Building the "GAlignmentPairs" object from the "first" and "last" mates
  # once the XS tag has been assigned to both

  r_XS <- GAlignmentPairs(r_first, r_last, strandMode = strandMode)

  return(r_XS)
}



