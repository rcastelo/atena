#' Sort by Query-name all BAM files in directory
#'
#' \code{sort_bam} allows to sort by Query-name (read ID) all binary 'BAM' files
#' located in the same directory. This function is based on the
#' \code{\link[Rsamtools]{sortBam}} function. This step is required
#' previous to the \code{\link{count_TE}} function when the suboptimal alignment
#' score tag is not specified in the 'BAM' files, in order to compute it.
#'
#' @param path_bam_files Character vector indicating the path in which the
#'   binary ‘BAM’ files to be sorted by the query-name (read ID) are located.
#' @param path_bam_files_sorted Character vector indicating the location where
#'   the resulting sorted BAM files will be created. This directory must be
#'   different from \code{path_bam_files}, otherwise the unsorted 'BAM' files
#'   will be overwritten.
#'
#' @return The function does not output any object. The resulting sorted files
#'   can be found in \code{path_bam_files_sorted}.
#'
#' @seealso
#' \code{\link{count_TE}},
#' \code{\link[Rsamtools]{sortBam}},
#'
#' @examples
#' \dontrun{path_bam_files <- "Documents/RNAseq/BAMfiles"
#' path_bam_files_sorted <- "Documents/RNAseq/BAMfiles_sorted"
#' sort_bam(path_bam_files, path_bam_files_sorted)}
#'
#' @export
#' @import Rsamtools



sort_bam <- function(path_bam_files,path_bam_files_sorted) {
  bam_files <- list.files(path_bam_files, pattern = "\\.bam$",
                          full.names = TRUE)
  bfl <- BamFileList(bam_files, asMates=TRUE)

  for (i in 1:length(bfl)) {
    file <- bfl[[i]]

    sortBam(file,
            destination = paste(path_bam_files_sorted,
                                sub(x = names(bfl[i]), pattern = ".bam",
                                    replacement = ""),
                                sep = "/"),
            byQname=TRUE)

  }
}
