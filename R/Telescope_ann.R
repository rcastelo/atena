#' Annotations of human transposable elements (HERVs and L1 LINEs) from Telescope
#' \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453}{(Bendall et al, 2019)}
#' 
#' \code{Telescope_ann()} provides the annotations of
#' transposable elements (human endogenous retroviruses (HERVs) and L1 long 
#' interspersed nuclear element (LINE)) created by Telescope authors
#' (\href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453}{Bendall et al, (2019)}). 
#' Transposable elements (TE) annotations in the GRCh38/hg38 human genome assembly include 72,169 ranges, 
#' corresponding to 28,513 TEs. For the GRCh37/hg19 assembly, there are 71,937
#' ranges corresponding to 28,314 TEs. The number difference between
#' ranges and TEs is due to the fact that some transcripts are disjoint, in 
#' order to exclude insertions of other element types. Telescope, however,
#' reports counts at the TE-level.
#' 
#' @param type (Default \code{all}). Character vector indicating the type of
#' elements to retrieve. Options are: \code{HERV}, which outputs annotations 
#' from HERVs, and \code{L1}, which gets the annotations of L1 elements. 
#' \code{all} option outputs the whole set of annotations, which includes
#' both HERVs and L1 elements.
#' 
#' @param version (Default \code{hg38}). Character string indicating the human
#' genome version of the annotations (\code{hg38} or \code{hg19}).
#' 
#' @return A \code{GRanges} object.
#' 
#' @source \url{https://github.com/mlbendall/telescope_annotation_db/blob/master/builds}
#'
#' @examples
#' Telescope_ann()
#'
#' @references
#' Bendall ML et al. Telescope: characterization of the retrotranscriptome by
#' accurate estimation of transposable element expression.
#' PLOS Computational Biology, 15:e1006453, 2019.
#' \url{https://doi.org/10.1371/journal.pcbi.1006453}
#'
#' @export
Telescope_ann <- function(type = c("all", "HERV", "L1"), version = c("hg38", "hg19")) {
  if (version == "hg38") {
    ann <- readRDS(file.path(system.file("extdata", package="atena"), "Telescope_ann_hg38.rds"))
  } else {
    ann <- readRDS(file.path(system.file("extdata", package="atena"), "Telescope_ann_hg19.rds"))
  }
  if (type == "HERV") {
    ann <- ann[mcols(ann)$source == "rmsk",]
  } else if (type == "L1") {
    ann <- ann[mcols(ann)$source == "l1base",]
  }
  ann
}
