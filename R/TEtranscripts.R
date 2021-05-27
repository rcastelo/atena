#' Build a TEtranscripts parameter object
#'
#' Build an object of the class \code{TEtranscriptsParam}
#'
#' @param bfl a character string vector of BAM filenames.
#'
#' @param teFeatures A \code{GRanges} or \code{GRangesList} object with the
#' TE annotated features to be quantified. Elements in this object should have
#' names, which will be used as a grouping factor for genomic ranges forming a
#' common locus, unless other metadata column names are specified in the
#' \code{aggregateby} parameter.
#'
#' @param aggregateby Character vector with column names in the annotation
#' to be used to aggregate quantifications. By default, this is an empty vector,
#' which means that the names of the input \code{GRanges} or \code{GRangesList}
#' object given in the \code{teFeatures} parameter will be used to aggregate
#' quantifications.
#'
#' @param geneFeatures A \code{GRanges} or \code{GRangesList} object with the
#' gene annotated features to be quantified. Following the TEtranscripts
#' algorithm, overlaps with unique reads are first tallied with respect to these
#' gene features.
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
#' will be considered to be overlapping an annotated feature as long as they
#' have a non-empty intersecting genomic ranges on the same strand, while when
#' \code{ignoreStrand = TRUE} the strand will not be considered.
#'
#' @param fragments (Default FALSE) A logical; applied to paired-end data only.
#' When \code{fragments=TRUE}, the read-counting method will also
#' count reads without mates, while when \code{fragments=FALSE} (default), those
#' reads will not be counted. For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}()}.
#' 
#' @param tolerance A positive numeric scalar storing the minimum tolerance
#' above which the SQUAREM algorithm (Du and Varadhan, 2020) keeps iterating.
#' Default is \code{1e-4} and this value is passed to the \code{tol} parameter
#' of the \code{\link[SQUAREM]{squarem}()} function.
#'
#' @param maxIter A positive integer scalar storing the maximum number of
#' iterations of the SQUAREM algorithm (Du and Varadhan, 2020). Default
#' is 100 and this value is passed to the \code{maxiter} parameter of the
#' \code{\link[SQUAREM]{squarem}()} function.
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
TEtranscriptsParam <- function(bfl, teFeatures, aggregateby=character(0),
                               geneFeatures=NA,
                               singleEnd=TRUE,
                               ignoreStrand=FALSE,
                               strandMode=1L,
                               fragments=FALSE,
                               tolerance=0.0001,
                               maxIter=100L) {

  bfl <- .checkBamFileListArgs(bfl, singleEnd, fragments)

  features <- .processFeatures(teFeatures, deparse(substitute(teFeatures)),
                               geneFeatures, deparse(substitute(geneFeatures)),
                               aggregateby)

  new("TEtranscriptsParam", bfl=bfl, features=features,
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
            cat(sprintf("# features (%s length %d): %s\n", class(object@features),
                        length(object@features),
                        ifelse(is.null(names(object@features)),
                               paste("on", .pprintnames(seqlevels(object@features))),
                               .pprintnames(names(object@features)))))
            cat(sprintf("# aggregated by: %s\n", ifelse(length(object@aggregateby) > 0,
                                                        paste(object@aggregateby, collapse=", "),
                                                        paste(class(object@features), "names"))))
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
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @aliases qtex
#' @aliases qtex,TEtranscriptsParam-method
#' @rdname qtex
setMethod("qtex", "TEtranscriptsParam",
          function(x, phenodata=NULL, mode=ovUnion, yieldSize=1e6L,
                   BPPARAM=SerialParam(progressbar=TRUE)) {
            .checkPhenodata(phenodata, length(x@bfl))

            cnt <- bplapply(x@bfl, .qtex_tetranscripts, ttpar=x, mode=mode,
                            yieldSize=yieldSize, BPPARAM=BPPARAM)
            cnt <- do.call("cbind", cnt)
            colData <- .createColumnData(cnt, phenodata)
            colnames(cnt) <- rownames(colData)

            features <- .consolidateFeatures(x, rownames(cnt))

            SummarizedExperiment(assays=list(counts=cnt),
                                 rowRanges=features,
                                 colData=colData)
          })

#' @importFrom stats setNames
#' @importFrom Rsamtools scanBamFlag ScanBamParam yieldSize yieldSize<-
#' @importFrom GenomicRanges width
#' @importFrom GenomicAlignments readGAlignments readGAlignmentsList
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom Matrix Matrix rowSums colSums t which
#' @importFrom SQUAREM squarem
#' @importFrom IRanges ranges
.qtex_tetranscripts <- function(bf, ttpar, mode, yieldSize=1e6L) {

  mode=match.fun(mode)
  readfun <- .getReadFunction(ttpar@singleEnd, ttpar@fragments)

  sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                         isDuplicate=FALSE,
                         isNotPassingQualityControls=FALSE)
  param <- ScanBamParam(flag=sbflags, tag="AS")

  if (is.null(ttpar@features$isTE)) {
    iste <- rep(TRUE, length(ttpar@features))
  } else {
    # gene <- ttpar@features[!ttpar@features$isTE & ttpar@features$type == "exon"]
    # names(gene) <- gene$gene_id
    # ttpar@features <- c(ttpar@features[ttpar@features$isTE], gene)
    iste <- as.vector(ttpar@features$isTE)
  }
  
  ov <- Hits(nLnode=0, nRnode=length(ttpar@features), sort.by.query=TRUE)
  alnreadids <- character(0)
  avgreadlen <- 0
  
  strand_arg <- "strandMode" %in% formalArgs(readfun)
  yieldSize(bf) <- yieldSize
  open(bf)
  while (length(alnreads <- do.call(readfun, c(list(file = bf), 
                                               list(param=param), 
                                               list(strandMode=ttpar@strandMode)[strand_arg], 
                                               list(use.names=TRUE))))) {
    avgreadlen <- avgreadlen + sum(width(ranges(alnreads)))
    alnreadids <- c(alnreadids, names(alnreads))
    thisov <- mode(alnreads, ttpar@features, ignoreStrand=ttpar@ignoreStrand)
    ov <- .appendHits(ov, thisov)
  }
  close(bf)

  ## getting the average fragment length
  nmappedreads <- length(alnreadids)
  avgreadlen <- avgreadlen / nmappedreads 

  ## count uniquely aligned-reads
  maskuniqaln <- !(duplicated(alnreadids) | duplicated(alnreadids, fromLast = TRUE))
  
  ## Adjusting quantification of multi-mapping reads overlapping common regions
  ## between genes and TEs: overlaps of multi-mapping reads mapping to genes 
  ## are removed
  if (!is.null(iste)) {
    ismulti <- !maskuniqaln[queryHits(ov)]
    iste_m <- iste[subjectHits(ov[ismulti])]
    iste_m <- split(x=iste_m, queryHits(ov[ismulti]))
    ovgenete <- unlist(lapply(iste_m, function(x) {length(unique(x)) != 1}))
    l <- unlist(lapply(iste_m, length))
    ovgenete <- rep(ovgenete, l)
    ov <- ov[-which(ismulti)[ovgenete & !unlist(iste_m)]]
  }
  
  ## fetch all different read identifiers from the overlapping alignments
  readids <- unique(alnreadids[queryHits(ov)])

  ## fetch all different transcripts from the overlapping alignments
  tx_idx <- sort(unique(subjectHits(ov)))

  ## build a matrix representation of the overlapping alignments
  ovalnmat <- .buildOvAlignmentsMatrix(ov, alnreadids, readids, tx_idx)


  # rsovalnmat <- rowSums(ovalnmat)
  # stopifnot(all(rsovalnmat > 0)) ## QC
  # maskuniqaln <- rsovalnmat == 1
  # rm(rsovalnmat)
  
  mt <- match(readids, alnreadids)
  
  if (!is.null(iste)) {
    ## Assigning unique reads mapping to an overlapping region between a TE and a 
    ## gene as gene counts
    ## which unique reads overlap both genes and TEs?
    istex <- as.vector(iste[tx_idx])
    idx <- (rowSums(ovalnmat[maskuniqaln[mt],istex]) > 0) & (rowSums(ovalnmat[maskuniqaln[mt],!istex]) > 0)
    ## Removing overlaps of unique reads with TEs if they also overlap a gene
    if (length(idx)>0) {
      ovalnmat[maskuniqaln[mt],][idx,istex] <- FALSE
    }
  }
  
  uniqcnt <- rep(0L, length(ttpar@features))
  uniqcnt[tx_idx] <- colSums(ovalnmat[maskuniqaln[mt], ])

  ## initialize vector of counts derived from multi-mapping reads
  cntvec <- rep(0, length(ttpar@features))

  if (sum(!maskuniqaln[mt]) > 0) { ## multi-mapping reads

    ## TEtranscripts doesn't use uniquely-aligned reads to inform
    ## the procedure of distributing multiple-mapping reads because,
    ## as explained in Jin et al. (2015) pg. 3594, "The unique-reads
    ## are not used as a prior for the initial abundance estimates in
    ## the EM procedure to reduce potential bias to certain TEs."
    ## for this reason, once counted, we discard unique alignments
    ovalnmat <- ovalnmat[!maskuniqaln[mt], ]
    readids <- readids[!maskuniqaln[mt]]
    ## the Qmat matrix stores row-wise the probability that read i maps to
    ## a transcript j, assume uniform probabilities by now
    Qmat <- Matrix(0, nrow=length(readids), ncol=length(tx_idx),
                   dimnames=list(readids, NULL))
    Qmat[which(ovalnmat, arr.ind=TRUE)] <- 1 ## Qmat is identical to ovalnmant except for rownames
    Qmat <- Qmat / rowSums(ovalnmat)

    ## Pi, corresponding to rho in Equations (1), (2) and (3) in
    ## Jin et al. (2015), stores probabilities of expression for each
    ## transcript, corrected for its effective length as defined
    ## in Eq. (1) of Jin et al. (2015)
    Pi <- colSums(Qmat)
    elen <- width(ttpar@features[tx_idx]) - avgreadlen + 1
    Pi <- .correctForTxEffectiveLength(Pi, elen)

    ## as specified in Jin et al. (2015), use the SQUAREM algorithm
    ## to achieve faster EM convergence
    emres <- squarem(par=Pi, Q=Qmat, elen=elen,
                     fixptfn=.ttFixedPointFun,
                     control=list(tol=ttpar@tolerance, maxiter=ttpar@maxIter))
    Pi <- emres$par
    Pi[Pi < 0] <- 0 ## Pi estimates are sometimes negatively close to zero
    Pi <- Pi / sum(Pi)

    ## use the estimated transcript expression probabilities
    ## to finally distribute ambiguously mapping reads
    probmassbyread <- as.vector(ovalnmat %*% Pi) 
    cntvecovtx <- rowSums(t(ovalnmat / probmassbyread) * Pi, na.rm=TRUE)
    cntvec[tx_idx] <- cntvecovtx
  }

  ## add multi-mapping and unique-mapping counts
  cntvec <- cntvec + uniqcnt
  names(cntvec) <- names(ttpar@features)

  ## aggregate TE quantifications if necessary
  if (length(ttpar@aggregateby) > 0) {
    f <- .factoraggregateby(ttpar@features[iste], ttpar@aggregateby)
    stopifnot(length(f) == length(cntvec[which(iste)])) ## QC
    cntvec_t <- tapply(cntvec[which(iste)], f, sum, na.rm=TRUE)
  }
  
  ## aggregating exon counts to genes
  if (!is.null(iste)) {
    fgene <- mcols(ttpar@features[!iste])[, "gene_id"]
    cntvec_g <- tapply(cntvec[which(!iste)], fgene, sum, na.rm=TRUE)
    cntvec <- c(cntvec_t, cntvec_g)
  } else {
    cntvec <- cntvec_t
  }

  ## TEtranscripts original implementation ultimately coerces fractional
  ## counts to integer 
  setNames(as.integer(cntvec), names(cntvec))
}

## private function .buildOvAlignmentsMatrix()
## builds a matrix representation of the overlapping alignments
## with reads on the rows and transcripts on the columns and
## a cell (i, j) set to TRUE if read i aligns to transcript j.
## therefore, when a read i aligns more than once to the same
## transcript j this is represented only once in the matrix
## parameters: ov - Hits object with the overlaps between
##                  alignments and features
##             arids - alignment read identifiers
##             rids - unique read identifiers
##             fidx - features index

#' @importFrom Matrix Matrix
#' @importFrom S4Vectors queryHits subjectHits
.buildOvAlignmentsMatrix <- function(ov, arids, rids, fidx) {
  oamat <- Matrix(FALSE, nrow=length(rids), ncol=length(fidx))
  mt1 <- match(arids[queryHits(ov)], rids)
  mt2 <- match(subjectHits(ov), fidx)
  oamat[cbind(mt1, mt2)] <- TRUE

  oamat
}

## private function .correctForTxEffectiveLength()
## corrects input transcript expression probabilities by
## the transcript effective length as defined in Eq. (1) of
## Jin et al. (2015)
## parameters: x - transcript expression probabilities
##             elen - effective length of each transcript
.correctForTxEffectiveLength <- function(x, elen) {
  x[elen > 0] <- x[elen > 0] / elen[elen > 0]
  ## x[elen <= 0] <- 0 ## (apparently this is done in the original Python code
  ## if (sum(x) > 0)   ##  but we don't do it here to avoid numerical instability)
    x <- x / sum(x)
  x
}

## private function .ttEstep()
## E-step of the EM algorithm of TEtranscripts
.ttEstep <- function(Q, Pi) {
  X <- t(t(Q) * Pi)
  X <- X / rowSums(X)
  X
}

## private function .ttEstep()
## M-step of the EM algorithm of TEtranscripts
.ttMstep <- function(X) {
  Pi <- colSums(X) / sum(X)
  Pi
}

## private function .ttFixedPointFun()
## fixed point function of the EM algorithm of TEtranscripts
.ttFixedPointFun <- function(Pi, Q, elen) {
  X <- .ttEstep(Q, Pi)
  Pi2 <- .ttMstep(X)
  Pi2 <- .correctForTxEffectiveLength(Pi2, elen)
  Pi2
}

## private function .ttLlh()
## log-likelihood function of the TEtranscripts model
## (currently not used)
.ttLlh <- function(X, Pi, Q) {
  z <- t(t(Q) * Pi)
  mask <- z == 0
  z[mask] <- NA
  sum(X * log(z), na.rm=TRUE)
}
