#' Build a Telescope parameter object
#'
#' Build an object of the class \code{TelescopeParam}.
#'
#' @param bfl A \code{BamFile} or \code{BamFileList} object, or a character
#' string vector of BAM filenames.
#'
#' @param teFeatures A \code{GRanges} or \code{GRangesList} object. Elements
#' in this object should have names, which will be used as a grouping factor
#' for ranges forming a common locus.(equivalent to "locus" column in 
#' Telescope), unless other metadata column names are specified in the
#' \code{aggregateby} parameter.
#' 
#' @param aggregateby Character vector with column names from the annotation
#' to be used to aggregate quantifications. By default, this is an empty vector,
#' which means that the names of the input \code{GRanges} or \code{GRangesList}
#' object given in the \code{teFeatures} parameter are used to aggregate
#' quantifications.
#' 
#' @param singleEnd (Default TRUE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @param strandMode (Default 1) Numeric vector which can take values 0, 1 or 2.
#' The strand mode is a per-object switch on
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' objects that controls the behavior of the strand getter. See
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' class for further detail. If \code{singleEnd = TRUE}, then \code{strandMode}
#' is ignored.
#' 
#' @param ignoreStrand (Default FALSE) A logical which defines if the strand
#' should be taken into consideration when computing the overlap between reads
#' and annotated features. When \code{ignoreStrand = FALSE}, an aligned read
#' is considered to be overlapping an annotated feature as long as they
#' have a non-empty intersecting genomic range on the same strand, while when
#' \code{ignoreStrand = TRUE} the strand is not considered.
#' 
#' @param fragments (Default FALSE) A logical; applied to paired-end data only.
#' When \code{fragments=FALSE} (default), the read-counting method only counts 
#' ‘mated pairs’ from opposite strands, while when \code{fragments=TRUE},
#' same-strand pairs, singletons, reads with unmapped pairs and other fragments 
#' are also counted. For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}()}.
#' 
#' @param pi_prior (Default 0) A positive integer scalar indicating the prior 
#' on pi. This is equivalent to adding n unique reads.
#'
#' @param theta_prior (Default 0) A positive integer scalar storing the prior 
#' on Q. Equivalent to adding n non-unique reads.
#'
#' @param em_epsilon (Default 1e-7) A numeric scalar indicating the EM 
#' Algorithm Epsilon cutoff.
#' 
#' @param maxIter A positive integer scalar storing the maximum number of
#' iterations of the EM SQUAREM algorithm (Du and Varadhan, 2020). Default
#' is 100 and this value is passed to the \code{maxiter} parameter of the
#' \code{\link[SQUAREM]{squarem}()} function.
#'
#'
#' @details
#' This is the constructor function for objects of the class
#' \code{TelescopeParam-class}. This type of object is the input to the
#' function \code{\link{qtex}()} for quantifying expression of transposable
#' elements, which will call the Telescope algorithm with this type of object.
#'
#' @return A \linkS4class{TelescopeParam} object.
#'
#' @examples
#' bamfiles <- list.files(system.file("extdata", package="atena"),
#'                        pattern="*.bam", full.names=TRUE)
#' annot <- Telescope_ann()
#' tspar <- TelescopeParam(bamfiles, annot)
#' tspar
#'
#' @references
#' Bendall et al. Telescope: characterization of the retrotranscriptome by
#' accurate estimation of transposable element expression.
#' PLOS Comp. Biol. 2019;15(9):e1006453. DOI:
#' \url{https://doi.org/10.1371/journal.pcbi.1006453}
#'
#' @importFrom methods is new
#' @importFrom Rsamtools BamFileList
#' @importFrom S4Vectors mcols
#' @export
TelescopeParam <- function(bfl, teFeatures, aggregateby=character(0),
                           singleEnd=TRUE, 
                           strandMode=1L, 
                           ignoreStrand=FALSE,
                           fragments=FALSE, 
                           pi_prior=0L, 
                           theta_prior=0L, 
                           em_epsilon=1e-7,
                           maxIter=100L) {
  bfl <- .checkBamFileListArgs(bfl, singleEnd, fragments)
  
  features <- .processFeatures(teFeatures, deparse(substitute(teFeatures)),
                               geneFeatures=NA, 
                               deparse(substitute(geneFeatures)),
                               aggregateby, 
                               aggregateexons = TRUE)
  
  new("TelescopeParam", bfl=bfl, features=teFeatures,
      aggregateby=aggregateby, singleEnd=singleEnd, ignoreStrand=ignoreStrand,
      strandMode=as.integer(strandMode), fragments=fragments,
      pi_prior=pi_prior, theta_prior=theta_prior, em_epsilon=em_epsilon,
      maxIter=as.integer(maxIter))
}

#' @param object A \linkS4class{TelescopeParam} object.
#'
#' @importFrom GenomeInfoDb seqlevels
#' @export
#' @aliases show,TelescopeParam-method
#' @rdname TelescopeParam-class
setMethod("show", "TelescopeParam",
          function(object) {
            cat(class(object), "object\n")
            cat(sprintf("# BAM files (%d): %s\n", length(object@bfl),
                        .pprintnames(names(object@bfl))))
            cat(sprintf("# features (%s length %d): %s\n", class(object@features),
                        length(object@features),
                        ifelse(is.null(names(object@features)),
                               paste("on", .pprintnames(seqlevels(object@features))),
                               .pprintnames(names(object@features)))))
            cat(sprintf("# aggregated by: %s\n", ifelse(length(object@aggregateby) > 0,
                                                        paste(object@aggregateby, collapse=", "),
                                                        paste(class(object@features), "names"))))
            cat(sprintf("# %s; %s",
                        ifelse(object@singleEnd, "single-end", "paired-end"),
                        ifelse(object@ignoreStrand, "unstranded", "stranded")))
            if (!object@ignoreStrand)
              cat(sprintf(" (strandMode=%d)", object@strandMode))
            if (!object@singleEnd)
              cat(sprintf("; %s",
                          ifelse(object@fragments, "counting properly paired, same-strand pairs, singletons, reads with unmapped pairs and other fragments",
                                 "counting properly paired reads")))
            cat("\n")
          })

#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @aliases qtex
#' @aliases qtex,TelescopeParam-method
#' @rdname qtex
setMethod("qtex", "TelescopeParam",
          function(x, phenodata=NULL, BPPARAM=SerialParam(progressbar=TRUE)) {
            .checkPhenodata(phenodata, length(x@bfl))

            cnt <- bplapply(x@bfl,
                            function(bf) {
                              tsexp <- basiliskRun(env=x@basiliskEnv,
                                                   fun=.qtex_telescope, bf=bf,
                                                   tspar=x)
                            }, BPPARAM=BPPARAM)
            cnt <- do.call("cbind", cnt)
            colData <- .createColumnData(cnt, phenodata)
            colnames(cnt) <- rownames(colData)

            SummarizedExperiment(assays=list(counts=cnt),
                                 rowRanges=x@features,
                                 colData=colData)
          })

#' @importFrom utils read.table
#' @importFrom reticulate import import_main
#' @importFrom BiocGenerics path
.qtex_telescope <- function(bf, tspar) {
  opts <- tspar@telescopeOptions
  opts$exp_tag <- gsub(".bam", "", basename(path(bf)))
  opts_str <- paste(paste0("--", names(opts)), sapply(opts, as.character))
  opts_str <- paste(gsub(" TRUE", "", opts_str), collapse=" ")
  opts_str_vec <- strsplit(opts_str, " ")[[1]]

  annfile <- file.path(tspar@telescopeOptions$outdir, "annotations.gtf")
  .exportTelescopeGTF(tspar@features, annfile)
  apmod <- reticulate::import("argparse")
  tsmod <- reticulate::import("telescope")
  pymain <- reticulate::import_main()
  pymain$sys$argv <- c("telescope", "assign", opts_str_vec, path(bf), annfile)
  desctxt <- "Tools for analysis of repetitive DNA elements"
  parser <- apmod$ArgumentParser(description=desctxt)
  VERSION <- tsmod[["_version"]]$VERSION
  parser$add_argument("--version", action="version",
                      version=VERSION, default=VERSION)
  subparsers <- parser$add_subparsers(help="sub-command help")
  assign_parser <- subparsers$add_parser("assign",
    description="Reassign ambiguous fragments that map to repetitive elements",
    formatter_class=apmod$ArgumentDefaultsHelpFormatter)
  tsmod$telescope_assign$IDOptions$add_arguments(assign_parser)
  assign_parser$set_defaults(func=tsmod$telescope_assign$run)
  pyargs <- parser$parse_args()
  pyargs$func(pyargs)
  dtf <- read.table(file.path(opts$outdir,
                              sprintf("%s-telescope_report.tsv", opts$exp_tag)),
                    header=TRUE)

  ## place quantifications in a vector matching the order of the annotations
  ## this also implies discarding the '__no_feature' quantification given by
  ## Telescope. note also that, with the exception of the '__no_feature',
  ## Telescope only outputs features that have a positive quantification.
  mt <- match(names(tspar@features), dtf$transcript)
  cntvec <- rep(0L, length=length(tspar@features))
  cntvec[!is.na(mt)] <- dtf$final_count[mt[!is.na(mt)]]
  names(cntvec) <- names(tspar@features)
  
  cntvec
}

## adapted from rtracklayer/R/gff.R
## private function to export a 'GRanges' or 'GRangesList' object to a GTF
## file in the GTF format expected by Telescope, which has the following two
## requirements:
##
## (1) one column must be called 'locus' and should contain
## values that group features forming part of the same locus
## (2) every line should end with a semicolon ';'
##
## arguments:
## gr : GRanges object
## fname : GTF filename
## src : value on the GTF source column

#' @importFrom methods is
#' @importFrom utils write.table
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols mcols<-
.exportTelescopeGTF <- function(gr, fname, src="Telescope") {

  if (is(gr, "GRangesList")) {
    gr <- unlist(gr)
  }

  con <- file(fname, open="wt")
  cat("", file=con) ## overwrite any existing file
  cat("##gff-version 2\n", file=con)

  seqs <- seqnames(gr)

  if (is.null(mcols(gr)$ID))
    mcols(gr)$ID <- names(gr)

  if (!is.null(mcols(gr)$source) && missing(src))
    src <- mcols(gr)$source
  else
    src <- rep(src, length(gr))

  features <- mcols(gr)$type
  if (is.null(features))
    features <- rep("feature", length(gr))

  scores <- mcols(gr)$score
  if (is.null(scores))
    scores <- rep(NA_real_, length(gr))

  strand <- strand(gr)
  if (is.null(strand))
    strand <- rep(strand(NA_character_), length(gr))
  strand[strand == "*"] <- NA_integer_

  phase <- mcols(gr)$phase
  if (is.null(phase)) {
    phase <- rep(NA_integer_, length(gr))
    if ("CDS" %in% features)
      warning("Annotation contains CDS features without phase information.")
  } else {
    if (anyNA(phase[features %in% "CDS"]))
      warning("Annotation contains CDS features with some phase information missing.")
  }

  if (is.null(mcols(gr)$locus))
      stop(".exportTelescopeGTF: input 'GRanges' object has no 'locus' metadata column.")

  tab <- data.frame(seqs, src, features, start(gr), end(gr),
                    scores, strand, phase)

  builtin <- c("type", "score", "phase", "source")
  custom <- setdiff(colnames(mcols(gr)), builtin)
  if (length(custom) > 0) {
    attrs <- mcols(gr)
    ## only deals with metadata columns that store atomic types, i.e.,
    ## not lists for instance.
    mdcols <- lapply(custom, function(name, attrs) {
                       x <- attrs[[name]]
                       res <- NULL
                       if (!is.atomic(x))
                         warning(sprintf("skipping non-atomic metadata column %s.", name))
                       else {
                         x_char <- as.character(x)
                         x_char <- sub(" *$", "", sub("^ *", "", as.character(x_char)))
                         x_char[is.na(x_char)] <- "\r"
                         if (!is.numeric(x))
                           x_char <- paste0("\"", x_char, "\"")
                         res <- paste(name, x_char, sep=" ")
                       }
                       res
                       }, attrs)
    mask <- sapply(lapply(mdcols,
                          function(x) if (is.logical(x) && !x) NULL else x),
                   is.null)
    mdcols[mask] <- NULL
    attrs <- as.data.frame(mdcols)
    attrs <- do.call(paste, c(attrs, sep="; "))
    attrs <- gsub("[^;]*?\r\"?(;|$)", "", attrs)
    attrs <- paste0(attrs, ";")
    attrs[nchar(attrs) == 0] <- NA
  }

  scipen <- getOption("scipen")
  options(scipen=100) ## prevent use of scientific notation
  on.exit(options(scipen=scipen))

  if (!is.null(attrs)) { # write out the rows with attributes first
    write.table(cbind(tab, attrs)[!is.na(attrs), ], con, sep="\t",
                na=".", quote=FALSE, col.names=FALSE, row.names=FALSE,
                append=TRUE)
    tab <- tab[is.na(attrs), ]
  }
  write.table(tab, con, sep="\t", na=".", quote=FALSE, col.names=FALSE,
              row.names=FALSE, append=TRUE)

  close(con)
}
