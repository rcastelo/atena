#' Build a Telescope parameter object
#'
#' Build an object of the class \code{TelescopeParam}.
#'
#' @param bfl A \code{BamFile} or \code{BamFileList} object, or a character
#' string vector of BAM filenames.
#'
#' @param annotations A \code{GRanges} object. Ranges in this object should
#' have names, which will be used as a grouping factor for ranges forming a
#' common locus.
#'
#' @param opts A \code{list} object specifying options to to pass to the
#' telescope algorithm, where a list element name should the option name and
#' its value should be the option value, which in the case of flags should be
#' be set as a logical \code{TRUE} or \code{FALSE} value. Defaults correspond
#' to those from the telescope software, with the exception of setting
#' \code{quiet=TRUE} by using \code{--quiet} in the call to telescope.
#' The following Telescope options cannot be set: \code{--attribute},
#' \code{--outdir}, \code{--exp_tag}, \code{--tempdir}, \code{--ncpu},
#' \code{--skip_em}. For a full documentation on Telescope options please
#' consult \url{https://github.com/mlbendall/telescope}.
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
#' \dontshow{basilisk.utils::installConda()}
#' bamfiles <- list.files(system.file("extdata", package="atena"),
#'                        pattern="*.bam", full.names=TRUE)
#' annot <- ERVmap_ann()
#' tspar <- TelescopeParam(bamfiles, annot)
#' tspar
#'
#' @references
#' Bendall et al. Telescope: characterization of the retrotranscriptome by
#' accurate estimation of transposable element expression.
#' PLOS Comp. Biol. 2019;15(9):e1006453. DOI:
#' \url{https://doi.org/10.1371/journal.pcbi.1006453}
#'
#' @importFrom basilisk BasiliskEnvironment basiliskStart basiliskStop
#' @export
TelescopeParam <- function(bfl, annotations, opts=list(quiet=TRUE)) {
  if (missing(bfl) || !class(bfl) %in% c("character", "BamFileList"))
    stop("argument 'bfl' should be either a string character vector of BAM file names or a 'BamFileList' object.")

  if (is.character(bfl)) {
    mask <- sapply(bfl, file.exists)
    if (any(!mask))
      stop(sprintf("the following input BAM files cannot be found:\n%s",
                   paste(paste("  ", bfl), collapse="\n")))
  }
  if (!is(bfl, "BamFileList"))
    bfl <- BamFileList(bfl)

  ## we may have to consider other types of objects than 'GRanges' to hold
  ## annotations (i.e., 'GRangesList', 'TxDb' ?)
  annotationsobjname <- deparse(substitute(annotations))
  env <- parent.frame()
  if (!exists(annotationsobjname))
    stop(sprintf("input GRranges object '%s' is not defined.", annotationsobjname))

  if (!is(annotations, "GRanges") && !is(annotations, "GRangesList"))
    stop(sprintf("annotations object '%s' should be either a 'GRanges' or a 'GRangesList' object.",
                 annotationsobjname))

  if (is.null(names(annotations)))
    stop(sprintf("the annotations object '%s' has no names.", annotationsobjname))

  gr <- annotations
  if (is(annotations, "GRangesList"))
    gr <- unlist(annotations)

  if ("locus" %in% colnames(mcols(gr)))
    warning(sprintf("the metadata column 'locus' in the annotations will be overwritten with the '%s' names.",
                    class(annotations)))

  gr$locus <- names(gr)

  if (is(annotations, "GRangesList"))
    annotations <- split(gr, names(gr))
  else ## is a GRanges object
    annotations <- gr

  if (!is.list(opts))
    stop("argument 'opts' should be a 'list' object.")

  if (length(opts) > 0) {
    bannedOpts <- c("attribute", "no_feature_key", "outdir", "exp_tag",
                    "tempdir", "ncpu", "skip_em")
    if (names(opts) %in% bannedOpts)
      stop(sprintf("The following Telescope options cannot be used from this package:\n  %s",
                   paste(bannedOpts, collapse=", ")))

    mask <- sapply(lapply(opts,
                          function(x) if (is.logical(x) && !x) NULL else x),
                   is.null)
    opts[mask] <- NULL
  }
  opts$outdir <- tempdir()

  pyenv <- BasiliskEnvironment("pyenv", pkgname="atena",
                                channels=c("conda-forge", "bioconda"),
                                packages=c("future==0.18.2", "pyyaml==5.3.1",
                                           "cython==0.29.7", "numpy==1.16.3",
                                           "scipy==1.2.1", "pysam>=0.16.0.1",
                                           "htslib==1.9", "intervaltree==3.0.2"),
                                path="telescope")
  cl <- basiliskStart(pyenv)
  tsversion <- basiliskRun(cl, function() {
                             tsmod <- reticulate::import("telescope")
                             tsmod[["_version"]]$VERSION
                           })
  basiliskStop(cl)

  new("TelescopeParam", bfl=bfl, annotations=annotations,
      basiliskEnv=pyenv, telescopeVersion=tsversion, telescopeOptions=opts)
}

#' @importFrom GenomeInfoDb seqlevels
#' @export
#' @aliases show,TelescopeParam-method
#' @rdname TelescopeParam-class
setMethod("show", "TelescopeParam",
          function(object) {
            cat(class(object), "object\n")
            cat(sprintf("# BAM files (%d): %s\n", length(object@bfl),
                        .pprintnames(names(object@bfl))))
            cat(sprintf("# annotations (%d): %s\n", length(object@annotations),
                        ifelse(is.null(names(object@annotations)),
                               paste("on", .pprintnames(seqlevels(object@annotations))),
                               .pprintnames(names(object@annotations)))))
            cat(sprintf("# Telescope version: %s\n", object@telescopeVersion))
            opts <- object@telescopeOptions
            opts_str <- paste(paste0("--", names(opts)), sapply(opts, as.character))
            opts_str <- paste(gsub(" TRUE", "", opts_str), collapse=" ")
            cat("# Telescope non-default options:\n")
            writeLines(paste0("#   ", strwrap(opts_str, width=60)))
          })

#' @param x A \code{TelescopeParam} object containing the configuration
#'        parameters to call the Telescope algorithm.
#'
#' @param phenodata A \code{data.frame} or \code{DataFrame} object storing
#'        phenotypic data to include in the resulting
#'        \code{SummarizedExperiment} object. If \code{phenodata} is set,
#'        its row names will become the column names of the resulting
#'        \linkS4class{SummarizedExperiment} object.
#'
#' @param BPPARAM An object of a \linkS4class{BiocParallelParam} subclass
#'        to configure the parallel execution of the code. By default,
#'        a \linkS4class{SerialParam} object is used, which does not use
#'        any parallelization, with the flag \code{progress=TRUE} to show
#'        progress through the calculations.
#'
#' @return A \linkS4class{SummarizedExperiment} object.
#'
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame
#' @export
#' @aliases qtex
#' @aliases qtex,TelescopeParam-method
#' @rdname qtex
setMethod("qtex", "TelescopeParam",
          function(x, phenodata=NULL, BPPARAM=SerialParam(progress=TRUE)) {
            if (!is.null(phenodata)) {
              if (nrow(phenodata) != length(x@bfl))
                stop("number of rows in 'phenodata' is different than the number of input BAM files in the input parameter object 'x'.")
              if (is.null(rownames(phenodata)))
                stop("'phenodata' has no row names.")
            }

            cnt <- bplapply(x@bfl,
                            function(bf) {
                              tsexp <- basiliskRun(env=x@basiliskEnv,
                                                   fun=.qtex_telescope, bf=bf,
                                                   tspar=x)
                            }, BPPARAM=BPPARAM)
            cnt <- do.call("cbind", cnt)

            colnames(cnt) <- gsub(".bam$", "", colnames(cnt))
            colData <- DataFrame(row.names=colnames(cnt))
            if (!is.null(phenodata)) {
              colData <- phenodata
              colnames(cnt) <- rownames(colData)
            }

            SummarizedExperiment(assays=list(counts=cnt),
                                 rowRanges=x@annotations,
                                 colData=colData)
          })

#' @importFrom reticulate import import_main
#' @importFrom BiocGenerics path
.qtex_telescope <- function(bf, tspar) {
  opts <- tspar@telescopeOptions
  opts$exp_tag <- gsub(".bam", "", basename(path(bf)))
  opts_str <- paste(paste0("--", names(opts)), sapply(opts, as.character))
  opts_str <- paste(gsub(" TRUE", "", opts_str), collapse=" ")
  opts_str_vec <- strsplit(opts_str, " ")[[1]]

  annfile <- file.path(tspar@telescopeOptions$outdir, "annotations.gtf")
  .exportTelescopeGTF(tspar@annotations, annfile)
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
  mt <- match(names(tspar@annotations), dtf$transcript)
  cntvec <- rep(0L, length=length(tspar@annotations))
  cntvec[!is.na(mt)] <- dtf$final_count[mt[!is.na(mt)]]
  names(cntvec) <- names(tspar@annotations)
  
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
  if (is.null(score))
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
  if (length(custom)) {
    attrs <- mcols(gr)
    attrs <- as.data.frame(sapply(custom, function(name) {
                                    x <- attrs[[name]]
                                    x_flat <- if (is(x, "List")) unlist(x, use.names=FALSE) else x
                                    x_char <- as.character(x_flat)
                                    x_char <- sub(" *$", "", sub("^ *", "", as.character(x_char)))
                                    if (is(x, "List")) {
                                      x_char[is.na(x_char)] <- "."
                                      x_char <- pasteCollapse(relist(x_char, x))
                                      x_char[lengths(x) == 0] <- NA
                                    }
                                    x_char[is.na(x_char)] <- "\r"
                                    if (!is.numeric(x_flat))
                                      x_char <- paste0("\"", x_char, "\"")
                                    paste(name, x_char, sep=" ")
                              }, simplify=FALSE))
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
