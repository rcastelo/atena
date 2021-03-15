#' Quantify transposable element expression
#'
#' The \code{qtex()} method quantifies transposable element expression.
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
#' @param phenodata A \code{data.frame} or \code{DataFrame} object storing
#'        phenotypic data to include in the resulting
#'        \code{SummarizedExperiment} object. If \code{phenodata} is set,
#'        its row names will become the column names of the resulting
#'        \linkS4class{SummarizedExperiment} object.
#'
#' @param mode One of the pre-defined overlapping methods such as
#'        \code{ovUnion()}, \code{ovIntersectionStrict} or a user-supplied
#'        overlapping function. For a user-supplied overlapping function, the
#'        input parameters must match those of the pre-defined methods and
#'        the function must return a \code{\link[S4Vectors]{Hits}} object with
#'        subject hits matching the annotated features. This parameter is
#'        analogous to the \code{mode} parameter of the
#'        \code{\link[GenomicAlignments]{summarizeOverlaps}()} function from
#'        the \code{GenomicAlignments} package.
#'
#' @param BPPARAM An object of a \linkS4class{BiocParallelParam} subclass
#'        to configure the parallel execution of the code. By default,
#'        a \linkS4class{SerialParam} object is used, which does not use
#'        any parallelization, with the flag \code{progress=TRUE} to show
#'        progress through the calculations.
#'
#' @return A \linkS4class{SummarizedExperiment} object.
#'
#' @details
#' Giving some \code{AtenaParam} object sub-class as input, the
#' \code{qtex()} method quantifies the expression of transposable
#' elements (TEs). The particular algorithm to perform the
#' quantification will be selected depending on the specific
#' sub-class of input \code{AtenaParam} object, see argument
#' \code{x} above.
#'
#' @seealso
#' \code{\link{ERVmapParam}}
#' \code{\link{TelescopeParam}}

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
#' @aliases qtex
#' @aliases qtex,AtenaParam-method
#' @rdname qtex
#' @name qtex
NULL

