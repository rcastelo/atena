#' Build a TEtranscripts parameter object
#'
#' Build an object of the class \code{TEtranscriptsParam}
#'
#' @param bfl A \code{BamFile} or \code{BamFileList} object, or a character
#' string vector of BAM filenames.
#'
#' @param annotations A \code{GRanges} or \code{GRangesList} object with the
#' annotations to be quantified. Elements in this object should have names,
#' which will be used as a grouping factor for genomic ranges forming a common
#' locus, unless other metadata column names are specified in the
#' \code{aggregateby} parameter.
#'
#' @param aggregateby Character vector with column names in the annotation
#' to be used to aggregate quantifications. By default, this is an empty vector,
#' which means that the names of the input \code{GRanges} or \code{GRangesList}
#' object given in the \code{annotations} parameter will be used to aggregate
#' quantifications.
#'
#' @param singleEnd (Default FALSE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @param strandMode (Default 1) Numeric vector which can take values 0, 1 or 2.
#'   The strand mode is a per-object switch on
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   objects that controls the behavior of the strand getter. See
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   class for further detail. If \code{singleEnd = TRUE}, then use either
#'   \code{strandMode = NULL} or do not specify the \code{strandMode} parameter.
#'
#' @param ignoreStrand (Default TRUE) A logical which defines if the strand
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
#' @param fragments (Default TRUE) A logical; applied to paired-end data only.
#' When \code{fragments=TRUE} (default), the read-counting method will also
#' count reads without mates, while when \code{fragments=FALSE} those reads
#' will not be counted. For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}()}.
#'
#' @details
#' This is the constructor function for objects of the class
#' \code{TEtranscriptsParam-class}. This type of object is the input to the
#' function \code{\link{qtex}()} for quantifying expression of transposable
#' elements using the TEtranscripts method
#' \href{https://doi.org/10.1093/bioinformatics/btv422}{Jin et al. (2015)}. The
#' TEtranscripts algorithm quantifies TE expression by using an EM algorithm
#' to optimally distribute ambiguously mapped reads.
#' 
#' @return A \linkS4class{TEtranscriptsParam} object.
#'
#' @examples
#' bamfiles <- list.files(system.file("extdata", package="atena"),
#'                        pattern="*.bam", full.names=TRUE)
#' annot <- ERVmap_ann()
#' ttpar <- TEtranscriptsParam(bamfiles, annot, singleEnd = TRUE)
#' ttpar
#'
#' @references
#' Jin Y et al. TEtranscripts: a package for including transposable elements
#' in differential expression analysis of RNA-seq datasets.
#' Bioinformatics. 2015;31(22):3593-3599. DOI:
#' \url{https://doi.org/10.1093/bioinformatics/btv422}
#'
#' @importFrom methods is new
#' @importFrom Rsamtools BamFileList
#' @export
TEtranscriptsParam <- function(bfl, annotations, aggregateby=character(0),
                               singleEnd=FALSE,
                               ignoreStrand=TRUE,
                               strandMode=1L,
                               fragments=TRUE,
                               tolerance=0.0001,
                               maxIter=100L) {
  if (missing(bfl) || !class(bfl) %in% c("character", "BamFileList"))
    stop("argument 'bfl' should be either a string character vector of BAM file names or a 'BamFileList' object")

  if (is.character(bfl)) {
    mask <- sapply(bfl, file.exists)
    if (any(!mask))
      stop(sprintf("The following input BAM files cannot be found:\n%s",
                   paste(paste("  ", bfl), collapse="\n")))
  }
  if (!is(bfl, "BamFileList"))
    bfl <- BamFileList(bfl, asMates=!singleEnd)

  annotationsobjname <- deparse(substitute(annotations))
  env <- parent.frame()
  if (!exists(annotationsobjname))
    stop(sprintf("input annotation object '%s' is not defined.", annotationsobjname))

  if (!is(annotations, "GRanges") && !is(annotations, "GRangesList"))
    stop(sprintf("annotations object '%s' should be either a 'GRanges' or a 'GRangesList' object.",
                 annotationsobjname))

  if (length(aggregateby) > 0)
    if (any(!aggregateby %in% colnames(mcols(annotations))))
        stop(sprintf("%s not in metadata columns of the annotations object.",
             paste(aggregateby[!aggregateby %in% colnames(mcols(annotations))])))

  if (is.null(names(annotations)) && length(aggregateby) == 0)
    stop(sprintf("the annotations object '%s' has no names and no aggregation metada columns have been specified.", annotationsobjname))

  if (is(annotations, "GRangesList"))
    annotations <- unlist(annotations)

  new("TEtranscriptsParam", bfl=bfl, annotations=annotations,
      aggregateby=aggregateby, singleEnd=singleEnd, ignoreStrand=ignoreStrand,
      strandMode=as.integer(strandMode), fragments=fragments,
      tolerance=tolerance, maxIter=as.integer(maxIter))
}

#' @param object A \linkS4class{TEtranscriptsParam} object.
#'
#' @importFrom GenomeInfoDb seqlevels
#' @export
#' @aliases show,TEtranscriptsParam-method
#' @rdname TEtranscriptsParam-class
setMethod("show", "TEtranscriptsParam",
          function(object) {
            cat(class(object), "object\n")
            cat(sprintf("# BAM files (%d): %s\n", length(object@bfl),
                        .pprintnames(names(object@bfl))))
            cat(sprintf("# annotations (%s length %d): %s\n", class(object@annotations),
                        length(object@annotations),
                        ifelse(is.null(names(object@annotations)),
                               paste("on", .pprintnames(seqlevels(object@annotations))),
                               .pprintnames(names(object@annotations)))))
            cat(sprintf("# aggregated by: %s\n", ifelse(length(object@aggregateby) > 0,
                                                        paste(object@aggregateby, collapse=", "),
                                                        paste(class(object@annotations), "names"))))
            cat(sprintf("# %s, %s",
                        ifelse(object@singleEnd, "single-end", "paired-end"),
                        ifelse(object@ignoreStrand, "unstranded", "stranded")))
            if (!object@ignoreStrand)
              cat(sprintf(" (strandMode=%d)", object@strandMode))
            if (!object@singleEnd)
              cat(sprintf(", %s",
                          ifelse(object@fragments, "counting each paired-end mate",
                                 "counting both paired-end mates")))
            cat("\n")
          })

#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRangesList
#' @export
#' @aliases qtex
#' @aliases qtex,TEtranscriptsParam-method
#' @rdname qtex
setMethod("qtex", "TEtranscriptsParam",
          function(x, phenodata=NULL, BPPARAM=SerialParam(progressbar=TRUE)) {
            if (!is.null(phenodata)) {
              if (nrow(phenodata) != length(x@bfl))
                stop("number of rows in 'phenodata' is different than the number of input BAM files in the input parameter object 'x'.")
              if (is.null(rownames(phenodata)))
                stop("'phenodata' has no row names.")
            }

            cnt <- bplapply(x@bfl, .qtex_tetranscripts, ttpar=x, BPPARAM=BPPARAM)

            cntmat <- do.call("cbind", cnt)
            colnames(cntmat) <- gsub(".bam$", "", names(x@bfl))
            colData <- DataFrame(row.names=colnames(cntmat))
            if (!is.null(phenodata)) {
              colData <- phenodata
              colnames(cntmat) <- rownames(colData)
            }

            annot <- x@annotations
            if (length(x@aggregateby) > 0) {
              f <- .factoraggregateby(x@annotations, x@aggregateby)
              annot <- GRangesList(split(annot, f))
            }

            SummarizedExperiment(assays=list(counts=cntmat),
                                 rowRanges=annot,
                                 colData=colData)
          })

#' @importFrom GenomicRanges mcols
f <- .factoraggregateby <- function(ann, aggby) {
  stopifnot(all(aggby %in% colnames(mcols(ann)))) ## QC
  spfstr <- paste(rep("%s", length(aggby)), collapse=":")
  f <- do.call("sprintf", c(spfstr, as.list(mcols(ann)[, aggby])))
  f
}

.correctForTxEffectiveLength <- function(x, elen) {
  x[elen > 0] <- x[elen > 0] / elen[elen > 0]
  ## x[elen <= 0] <- 0 ## (apparently this is done in the Python code
  x <- x / sum(x)      ##  but we don't do it here to avoid numerical instability)
  x
}

.ttEstep <- function(Q, Pi) {
  X <- t(t(Q) * Pi)
  X <- X / rowSums(X)
  X
}

.ttMstep <- function(X) {
  Pi <- colSums(X) / sum(X)
  Pi
}

.ttFixedPointFun <- function(Pi, Q, elen) {
  X <- .ttEstep(Q, Pi)
  Pi2 <- .ttMstep(X)
  Pi2 <- .correctForTxEffectiveLength(Pi2, elen)
  Pi2
}

.llh <- function(X, Pi, Q) {
  z <- t(t(Q) * Pi)
  mask <- z == 0
  z[mask] <- NA
  sum(X * log(z), na.rm=TRUE)
}

#' @importFrom Rsamtools scanBamFlag ScanBamParam yieldSize asMates
#' @importFrom GenomicRanges width
#' @importFrom GenomicAlignments readGAlignments readGAlignmentsList qwidth
#' @importFrom GenomicAlignments readGAlignmentPairs first last
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix Matrix rowSums colSums t which
#' @importFrom SQUAREM squarem
.qtex_tetranscripts <- function(bf, ttpar) {
  sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                         isDuplicate=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isSecondary=TRUE)

  param <- ScanBamParam(flag=sbflags, tag="AS")

  open(bf)
  ## use.names=TRUE is important to get the read names
  multialignments <- NULL
  if (ttpar@singleEnd)
    multialignments <- readGAlignments(bf, use.names=TRUE, param=param)
  else {
    if (ttpar@fragments)
      multialignments <- readGAlignmentsList(bf, use.names=TRUE, param=param)
    else
      multialignments <- readGAlignmentPairs(bf, use.names=TRUE, param=param)
  }
  close(bf)

  teann <- ttpar@annotations

  ## find overlaps between multimapping reads and annotations
  ovmulti <- findOverlaps(multialignments, teann,
                          ignore.strand=ttpar@ignoreStrand)

  ## fetch all different read names from the overlapping alignments
  read_names <- unique(names(multialignments)[queryHits(ovmulti)])

  ## fetch all different transcripts from the overlapping alignments
  tx_idx <- unique(subjectHits(ovmulti))

  ## build a matrix representation of the overlapping alignments
  ## with reads on the rows and transcripts on the columns and
  ## a cell (i, j) set to TRUE if read i aligns to transcript j
  ## when a read i aligns more than once to the same transcript j
  ## this is represented only once in the matrix
  ovalnmat <- Matrix(FALSE, nrow=length(read_names), ncol=length(tx_idx),
                                        dimnames=list(read_names, NULL))
  mt1 <- match(names(multialignments)[queryHits(ovmulti)], read_names)
  mt2 <- match(subjectHits(ovmulti), tx_idx)
  ovalnmat[cbind(mt1, mt2)] <- TRUE

  ## the Qmat matrix stores row-wise the probability that read i maps to
  ## a transcript j, assume uniform probabilities by now
  Qmat <- Matrix(0, nrow=length(read_names), ncol=length(tx_idx),
                 dimnames=list(read_names, NULL))
  Qmat[which(ovalnmat, arr.ind=TRUE)] <- 1
  Qmat <- Qmat / rowSums(ovalnmat)

  ## figure out average overlapping-read length
  avgReadLength <- NA_integer_
  if (length(read_names) <= 1000)
    avgReadLength <- mean(width(ranges(multialignments[read_names])))
  else {
    sam <- sample(1:length(read_names), size=1000, replace=TRUE)
    avgReadLength <- mean(width(ranges(multialignments[read_names[sam]])))
  }

  ## Pi, corresponding to rho in Equations (1), (2) and (3) in
  ## Jin et al. (2015), stores probabilities of expression for each
  ## transcript, corrected for its effective length as defined
  ## in Eq. (1) of Jin et al. (2015)
  Pi <- colSums(Qmat)
  elen <- width(teann[tx_idx]) - avgReadLength + 1
  Pi <- .correctForTxEffectiveLength(Pi, elen)

  ## as specified in Jin et al. (2015), use the SQUAREM algorithm
  ## to achieve faster EM convergence
  emres <- squarem(p=Pi, Q=Qmat, elen=elen,
                   fixptfn=.ttFixedPointFun,
                   control=list(tol=ttpar@tolerance, maxiter=ttpar@maxIter))
  Pi <- emres$par
  Pi[Pi < 0] <- 0 ## Pi estimates are sometimes negatively close to zero
  Pi <- Pi / sum(Pi)

  ## use the estimated transcript expression probabilities
  ## to finally distribute ambiguously mapping reads
  probmassbyread <- as.vector(ovalnmat %*% Pi)
  cntvecovtx <- rowSums(t(ovalnmat / probmassbyread) * Pi, na.rm=TRUE)
  cntvec <- rep(0, length(teann))
  cntvec[tx_idx] <- cntvecovtx

  ## add unique aligned reads discarding those overlapping more than
  ## one feature

  sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                         isDuplicate=FALSE,
                         isNotPassingQualityControls=FALSE,
                         isSecondary=FALSE)

  param <- ScanBamParam(flag=sbflags, tag="AS")

  open(bf)
  ## use.names=TRUE is important to get the read names
  uniqalignments <- NULL
  if (ttpar@singleEnd)
    uniqalignments <- readGAlignments(bf, use.names=TRUE, param=param)
  else {
    if (ttpar@fragments)
      uniqalignments <- readGAlignmentsList(bf, use.names=TRUE, param=param)
    else
      uniqalignments <- readGAlignmentPairs(bf, use.names=TRUE, param=param)
  }
  close(bf)

  ## find overlaps of unique-mapping alignments with TE annotations
  ovuniq <- findOverlaps(uniqalignments, teann, ignore.strand=TRUE)

  ## discard unique-mapping alignments overlapping more than one TE
  whuniqalnTEs <- which(countSubjectHits(ovuniq) == 1L)
  ovuniq <- ovuniq[subjectHits(ovuniq) %in% whuniqalnTEs]
  uniqcntvec <- countSubjectHits(ovuniq)

  ## add multi-mapping and unique-mapping counts
  cntvec <- cntvec + uniqcntvec

  ## aggregate quantifications if necessary
  if (length(ttpar@aggregateby) > 0) {
    f <- .factoraggregateby(teann, ttpar@aggregateby)
    stopifnot(length(f) == length(cntvec)) ## QC
    cntvec <- tapply(cntvec, f, sum, na.rm=TRUE)
  }

  as.integer(cntvec)
}
