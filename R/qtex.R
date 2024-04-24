#' Quantify transposable element expression
#'
#' The \code{qtex()} method quantifies transposable element expression.
#' 
#' @param x An \code{QuantifyParam} object of one of the following
#' subclasses:
#' \itemize{
#'   \item A \code{TEtranscriptsParam} object built using the constructor
#'         function \code{\link{TEtranscriptsParam}()}. This object will
#'         trigger \code{qtex()} to use the quantification algorithm by
#'         Jin et al. (2015).
#'   \item A \code{ERVmapParam} object built using the constructor
#'         function \code{\link{ERVmapParam}()}. This object will
#'         trigger \code{qtex()} to use the quantification algorithm by
#'         Tokuyama et al. (2018).
#'   \item A \code{TelescopeParam} object built using the constructor
#'         function \code{\link{TelescopeParam}()}. This object will
#'         trigger \code{qtex()} to use the quantification algorithm by
#'         Bendall et al. (2019).
#'   \item An \code{atenaParam} object built using the constructor
#'         function \code{\link{atenaParam}()}. This object will
#'         trigger \code{qtex()} to use a quantification algorithm
#'         specifically developed in this package.
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
#' @param yieldSize Field inherited from \code{\link[Rsamtools]{BamFile}}.
#'        The method for signature \code{\link{ERVmapParam}()} reads the BAM
#'        file by chunks. \code{yieldSize} represents the number of records
#'        (chunk size) to yield each time the file is read.
#'
#' @param auxiliaryFeatures (Default \code{FALSE}). It only applies when `x` is
#'        a [`TelescopeParam`] or an [`atenaParam`] object. When \code{TRUE},
#'        auxiliary features created during expression quantification are also
#'        returned in the [`SummarizedExperiment`] object.
#'
#' @param verbose (Default 1). When \code{verbose} > 1, detailed information on
#'        the quantification steps is provided. Warnings are always present
#'        regardless of the value of \code{verbose}.
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
#' \code{\link{TEtranscriptsParam}}
#' \code{\link{ERVmapParam}}
#' \code{\link{TelescopeParam}}
#'
#' @references
#' Jin Y et al. TEtranscripts: a package for including transposable elements
#' in differential expression analysis of RNA-seq datasets.
#' Bioinformatics. 2015;31(22):3593-3599. DOI:
#' \url{https://doi.org/10.1093/bioinformatics/btv422}
#'
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
#' rmskLTR <- getLTRs(rmskat, relLength=0.8,
#'                    fullLength=TRUE,
#'                    partial=TRUE,
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
#' ## quantify expression
#' qts <- qtex(tspar)
NULL
