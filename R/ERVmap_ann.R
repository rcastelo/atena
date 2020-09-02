#' Annotations of human endogenous retroviruses from ERVmap in
#' \href{https://doi.org/10.1073/pnas.1814589115}{Tokuyama et al. (2018)}
#'
#' Annotations of human endogenous retroviruses (HERVs) published in
#' \href{https://doi.org/10.1073/pnas.1814589115}{Tokuyama et al. (2018)} for
#' the ERVmap pipeline. These include 3220 curated autonomous HERVs that mirror
#' full-length proviruses.
#'
#' @format A \code{GRanges} object with 3220 ranges and 2 metadata columns:
#' \describe{
#'   \item{seqnames}{chromosome in which the HERV is located}
#'   \item{ranges}{start and end positions of the HERV in the chromosome}
#'   \item{strand}{strand of the HERV}
#'   \item{name}{name of the HERV. This is the published name except for, the
#'   HERVs without a name in the original studies, in which a numerical ID
#'   following the family name of the HERV is assigned as the name}
#'   \item{score}{score for the HERV annotation}
#' }
#' @source \url{https://github.com/mtokuyama/ERVmap/blob/master/ERVmap.bed}
#'
#' @references
#' Tokuyama M et al. ERVmap analysis reveals genome-wide transcription of human
#' endogenous retroviruses. PNAS. 2018;115(50):12565-12572. DOI:
#' \url{https://doi.org/10.1073/pnas.1814589115}
#'
"ERVmap_ann"
