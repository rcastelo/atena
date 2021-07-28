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
#' @param geneFeatures A \code{GRanges} or \code{GRangesList} object with the
#' gene annotated features to be quantified. The TEtranscripts approach for
#' gene expression quantification is used, in which overlaps with unique reads 
#' are first tallied with respect to these gene features whereas multi-mapping
#' reads are preferentially assigned to TEs. Elements should have names 
#' indicating the gene name/id. In case that \code{geneFeatures} contains a 
#' metadata column named \code{type}, only the elements with 
#' \code{type} = \code{exon} are considered for the analysis. Then, exon counts 
#' are summarized to the gene level.
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
#' elements, which will call the Telescope algorithm 
#' \href{https://doi.org/10.1371/journal.pcbi.1006453}{Bendall et al. (2019)} 
#' with this type of object.
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
                           geneFeatures=NA,
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
                               geneFeatures, deparse(substitute(geneFeatures)),
                               aggregateby, aggregateexons = TRUE)
  
  new("TelescopeParam", bfl=bfl, features=features,
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

#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @aliases qtex
#' @aliases qtex,TelescopeParam-method
#' @rdname qtex
setMethod("qtex", "TelescopeParam",
          function(x, phenodata=NULL, mode=ovUnion, yieldSize=1e6L,
                   BPPARAM=SerialParam(progressbar=TRUE)) {
            .checkPhenodata(phenodata, length(x@bfl))
            
            cnt <- bplapply(x@bfl, .qtex_telescope, tspar=x, mode=mode,
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
#' @importFrom BiocGenerics path
#' @importFrom GenomicAlignments readGAlignments readGAlignmentsList
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom Matrix Matrix rowSums colSums t which
#' @importFrom SQUAREM squarem
.qtex_telescope <- function(bf, tspar, mode, yieldSize=1e6L) {
  mode=match.fun(mode)
  
  readfun <- .getReadFunction(tspar@singleEnd, tspar@fragments)
  sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                         isDuplicate=FALSE,
                         isNotPassingQualityControls=FALSE)
  param <- ScanBamParam(flag=sbflags, tag="AS")
  
  iste <- as.vector(attributes(tspar@features)$isTE[,1])
  
  if (any(duplicated(names(tspar@features[iste])))) {
    stop(".qtex_telescope: transposable element annotations do not contain unique names for each element")
  }
  
  ov <- Hits(nLnode=0, nRnode=length(tspar@features), sort.by.query=TRUE)
  alnreadids <- character(0)
  avgreadlen <- integer()
  asvalues <- integer()
  
  strand_arg <- "strandMode" %in% formalArgs(readfun)
  yieldSize(bf) <- yieldSize
  open(bf)
  while (length(alnreads <- do.call(readfun, c(list(file = bf), 
                                               list(param=param), 
                                               list(strandMode=tspar@strandMode)[strand_arg], 
                                               list(use.names=TRUE))))) {
    avgreadlen <- c(avgreadlen, width(ranges(alnreads)))
    alnreadids <- c(alnreadids, names(alnreads))
    asvalues <- c(asvalues, mcols(alnreads)$AS)
    thisov <- suppressWarnings(mode(alnreads, tspar@features, 
                                    ignoreStrand=tspar@ignoreStrand))
    ov <- .appendHits(ov, thisov)
  }
  close(bf)
  
  ## create a logical mask for uniquely aligned-reads
  maskuniqaln <- !(duplicated(alnreadids) | 
                     duplicated(alnreadids, fromLast = TRUE))
  
  ## fetch all different read identifiers from the overlapping alignments
  readids <- unique(alnreadids[queryHits(ov)])
  
  ## fetch all different transcripts from the overlapping alignments
  tx_idx <- sort(unique(subjectHits(ov)))
  istex <- as.vector(iste[tx_idx])
  
  ## build a matrix representation of the overlapping alignments
  ovalnmat <- .buildOvAlignmentsMatrix(ov, alnreadids, readids, tx_idx)
  
  ## store logical mask for multi-mapping reads for only those present in the 
  ## ov matrix
  mt <- match(readids, alnreadids)
  maskmulti <- !maskuniqaln[mt]
  
  if (!all(iste)) {
    ## Correcting for preference of unique/multi-mapping reads to genes or TEs, respectively
    ovalnmat <- .correctPreference(ovalnmat, maskuniqaln, mt, istex)
    stopifnot(!any(rowSums(ovalnmat[,istex]) > 0 & rowSums(ovalnmat[,!istex]) > 0))
  }
  
  ## initialize vector of counts derived from multi-mapping reads
  cntvec <- rep(0L, length(tspar@features))
  
  # Getting counts from reads overlapping only one feature (this excludes
  # unique reads mapping to two or more overlapping features)
  ovunique <- rowSums(ovalnmat) == 1
  cntvec[tx_idx] <- colSums(ovalnmat[ovunique,])
  
  # --- EM-step --- 
  alnreadidx <- match(alnreadids, readids)
  rd_idx <- sort(unique(alnreadidx[queryHits(ov)]))
  asvalues <- (asvalues - min(asvalues)) / (max(asvalues) - min(asvalues))
  QmatTS <- .buildOvValuesMatrixTS(ov, asvalues, alnreadidx, rd_idx, tx_idx)
  QmatTS <- QmatTS / rowSums(QmatTS)
  QmatTS[is.na(rowSums(QmatTS)),] <- rep(0, )
  # Telescope [Bendall et al. (2019)) defines the initial π estimate uniformly.
  PiTS <- rep(1 / length(tx_idx), length(tx_idx))
  
  # Telescope also defines an additional reassignment parameter (θ) as uniform
  Theta <- rep(1 / length(tx_idx), length(tx_idx))
  
  ## The SQUAREM algorithm to run the EM procedure
  a <- tspar@pi_prior # 0
  b <- tspar@theta_prior # 0
  GlobalTheta <<- Theta
  tsres <- squarem(par=PiTS, Q=QmatTS, maskmulti=maskmulti, a=a, b=b,
                   fixptfn=.tsFixedPointFun,
                   control=list(tol=tspar@em_epsilon, maxiter=tspar@maxIter))
  # maskmulti_ovmulti <- maskmulti | !ovunique
  # tsres <- squarem(par=PiTS, Q=QmatTS, maskmulti=maskmulti_ovmulti, a=a, b=b,
  #                  fixptfn=.tsFixedPointFun,
  #                  control=list(tol=tspar@em_epsilon, maxiter=tspar@maxIter))
  PiTS <- tsres$par
  PiTS[PiTS < 0] <- 0 ## Pi estimates are sometimes negatively close to zero
  # --- end EM-step ---
  
  # Implementation for reassign_mode=exclude
  X <- .tsEstep(QmatTS, GlobalTheta, maskmulti, PiTS)
  maxbyrow <- rowMaxs(X)
  Xind <- X == maxbyrow
  nmaxbyrow <- rowSums(Xind)
  #cntvec[tx_idx][istex] <- colSums(Xind[nmaxbyrow == 1, ])
  cntvec[tx_idx] <- colSums(Xind[nmaxbyrow == 1 & !ovunique, ]) # ovunique reads have already been count. Do not differenciate between TEs and genes with istex because in Telescope genes are also included in the EMstep.
  
  ## aggregate TE quantifications if necessary
  if (length(tspar@aggregateby) > 0) {
    f <- .factoraggregateby(tspar@features[iste], tspar@aggregateby)
    stopifnot(length(f) == length(cntvec)) ## QC
    cntvec <- tapply(cntvec, f, sum, na.rm=TRUE)
  }
  
  setNames(as.integer(cntvec), names(cntvec))
  
}



## private function .tsEstep()
## E-step of the EM algorithm of Telescope
.tsEstep <- function(Q, Theta, maskmulti, Pi) {
  X <- t(t(Q) * Pi)
  X[maskmulti, ] <- t(t(X[maskmulti, ]) * Theta)
  X <- X[rowSums(X)>0,, drop=FALSE]
  X <- X / rowSums(X)
  X
}

## private function .tsMstepPi()
## M-step of the EM algorithm of Telescope
.tsMstepPi <- function(X, a) {
  Pi <- (colSums(X) + a) / (sum(X)+a*ncol(X))
  Pi
}

## private function .tsMstepTheta()
## Update the estimate of the MAP value of θ
.tsMstepTheta <- function(X, maskmulti, b) {
  Theta <- (colSums(X[maskmulti, , drop=FALSE]) + b) / (sum(maskmulti) + b*ncol(X))
}

## private function .tsFixedPointFun()
## fixed point function of the EM algorithm of Telescope
.tsFixedPointFun <- function(Pi, Q, maskmulti, a, b) {
  Theta <- GlobalTheta
  X <- .tsEstep(Q, Theta, maskmulti, Pi)
  Pi2 <- .tsMstepPi(X, a)
  Theta2 <- .tsMstepTheta (X, maskmulti, b)
  GlobalTheta <<- Theta2 ## need to work with a global variable
  
  Pi2
}

## private function .buildOvValuesMatrixTS()
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix Matrix
.buildOvValuesMatrixTS <- function(ov, values, aridx, ridx, fidx) {
  ovmat <- Matrix(do.call(class(values), list(1)),
                  nrow=length(ridx), ncol=length(fidx))
  mt1 <- match(aridx[queryHits(ov)], ridx)
  mt2 <- match(subjectHits(ov), fidx)
  ovmat[cbind(mt1, mt2)] <- values[queryHits(ov)]
  ovmat
}





