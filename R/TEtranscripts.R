#' Build a TEtranscripts parameter object
#'
#' Build an object of the class \code{TEtranscriptsParam}
#'
#' @param bfl a character string vector of BAM file names.
#'
#' @param teFeatures A \code{GRanges} or \code{GRangesList} object with the
#' TE annotated features to be quantified. Elements in this object should have
#' names, which are used as a grouping factor for genomic ranges forming a
#' common locus, unless other metadata column names are specified in the
#' \code{aggregateby} parameter.
#'
#' @param aggregateby Character vector with column names from the annotation
#' to be used to aggregate quantifications. By default, this is an empty
#' vector, which means that the names of the input \code{GRanges} or
#' \code{GRangesList} object given in the \code{teFeatures} parameter are used
#' to aggregate quantifications.
#' 
#' @param ovMode Character vector indicating the overlapping mode. Available
#' options are: "ovUnion" (default) and "ovIntersectionStrict",
#' which implement the corresponding methods from HTSeq
#' (\url{https://htseq.readthedocs.io/en/release_0.11.1/count.html}).
#' Ambiguous alignments (alignments overlapping > 1 feature) are addressed
#' as in the original TEtranscripts method.
#'
#' @param geneFeatures A \code{GRanges} or \code{GRangesList} object with the
#' gene annotated features to be quantified. Following the TEtranscripts
#' algorithm, overlaps with unique reads are first tallied with respect to
#' these gene features. Elements should have names indicating the gene name/id.
#' In case that \code{geneFeatures} contains a metadata column named
#' \code{type}, only the elements with \code{type} = \code{exon} are considered
#' for the analysis. Then, exon counts are summarized to the gene level.
#'
#' @param singleEnd (Default TRUE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @param strandMode (Default 1) Numeric vector which can take values 0, 1 or
#' 2.
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
#' @param fragments (Default TRUE) A logical; applied to paired-end data only.
#' In both cases (\code{fragments=FALSE} and \code{fragments=TRUE}), the
#' read-counting method discards not properly paired reads. Moreover,
#' when \code{fragments=FALSE}, only non-ambiguous properly paired reads are
#' counted. When \code{fragments=TRUE}, ambiguous reads are also counted 
#' (see "Pairing criteria" in \code{\link[GenomicAlignments]{readGAlignments}()}). 
#' \code{fragments=TRUE} is equivalent to
#' the behavior of the TEtranscripts algorithm. For further details see
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
#' TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
#'                     package="atena"))
#' ttpar <- TEtranscriptsParam(bamfiles, teFeatures = TE_annot,
#'                             singleEnd = TRUE, ignoreStrand=TRUE,
#'                             aggregateby = c("repName"))
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
#' @rdname TEtranscriptsParam-class
TEtranscriptsParam <- function(bfl, teFeatures, aggregateby=character(0),
                                ovMode="ovUnion",
                                geneFeatures=NA,
                                singleEnd=TRUE,
                                ignoreStrand=FALSE,
                                strandMode=1L,
                                fragments=TRUE,
                                tolerance=0.0001,
                                maxIter=100L) {

    bfl <- .checkBamFileListArgs(bfl, singleEnd, fragments)

    if (!ovMode %in% c("ovUnion","ovIntersectionStrict"))
      stop("'ovMode' should be one of 'ovUnion', 'ovIntersectionStrict'")
    
    features <- .processFeatures(teFeatures, deparse(substitute(teFeatures)),
                                geneFeatures, deparse(substitute(geneFeatures)),
                                aggregateby, aggregateexons = FALSE)

    new("TEtranscriptsParam", bfl=bfl, features=features,
        aggregateby=aggregateby, ovMode=ovMode,
        singleEnd=singleEnd, ignoreStrand=ignoreStrand,
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
            cat(sprintf("# features (%s length %d): %s\n",
                        class(object@features),
                        length(object@features),
                        ifelse(is.null(names(object@features)),
                                paste("on",
                                    .pprintnames(seqlevels(object@features))),
                                .pprintnames(names(object@features)))))
            cat(sprintf("# aggregated by: %s\n",
                        ifelse(length(object@aggregateby) > 0,
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
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @aliases qtex
#' @aliases qtex,TEtranscriptsParam-method
#' @rdname qtex
setMethod("qtex", "TEtranscriptsParam",
        function(x, phenodata=NULL, mode=ovUnion, yieldSize=1e6L,
                    BPPARAM=SerialParam(progressbar=TRUE)) {
            .checkPhenodata(phenodata, length(x@bfl))

            cnt <- bplapply(x@bfl, .qtex_tetranscripts, ttpar=x, mode=x@ovMode,
                            yieldSize=yieldSize, BPPARAM=BPPARAM)
            cnt <- do.call("cbind", cnt)
            colData <- .createColumnData(cnt, phenodata)
            colnames(cnt) <- rownames(colData)

            features <- .consolidateFeatures(x, rownames(cnt))

            SummarizedExperiment(assays=list(counts=cnt),
                                rowRanges=features,
                                colData=colData)
        })

#' @importFrom stats setNames aggregate median
#' @importFrom Rsamtools scanBamFlag ScanBamParam yieldSize yieldSize<-
#' @importFrom GenomicRanges width
#' @importFrom GenomicAlignments readGAlignments readGAlignmentsList
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom Matrix Matrix rowSums colSums t which
#' @importFrom IRanges ranges
.qtex_tetranscripts <- function(bf, ttpar, mode, yieldSize=1e6L) {
    mode <- match.fun(mode)
    readfun <- .getReadFunction(ttpar@singleEnd, ttpar@fragments)
    sbflags <- scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE,
                           isNotPassingQualityControls=FALSE,
                           isProperPair=TRUE)
    param <- ScanBamParam(flag=sbflags, what="flag", tag="AS")
    
    iste <- as.vector(attributes(ttpar@features)$isTE[,1])
    if (any(duplicated(names(ttpar@features[iste])))) {
        stop(".qtex_tetranscripts: transposable element annotations do not contain unique names for each element")
    }
    ov <- Hits(nLnode=0, nRnode=length(ttpar@features), sort.by.query=TRUE)
    alnreadids <- character(0)
    avgreadlen <- integer(0)
    d <- numeric(0)
    salnmask <- logical(0)
    strand_arg <- "strandMode" %in% formalArgs(readfun)
    yieldSize(bf) <- yieldSize
    open(bf)
    while (length(alnreads <- do.call(readfun, 
                                c(list(file = bf), list(param=param),
                                list(strandMode=ttpar@strandMode)[strand_arg],
                                list(use.names=TRUE))))) {
        avgreadlen <- c(avgreadlen, .getAveLen(ttpar, alnreads))
        if (is(alnreads, "GAlignmentPairs")) {
            d <- c(d, abs(start(first(alnreads)) - start(second(alnreads))))
        } else if (is(alnreads, "GAlignmentsList")) {
            d <- c(d, max(abs(diff(start(alnreads)))))
            d[d < 0] <- 0
        }
        salnmask <- c(salnmask, any(.secondaryAlignmentMask(alnreads)))
        alnreadids <- c(alnreadids, names(alnreads))
        thisov <- mode(alnreads, ttpar@features,
                        ignoreStrand=ttpar@ignoreStrand, inter.feature=FALSE)
        ov <- .appendHits(ov, thisov)
    }
    # close(bf)
    on.exit(close(bf))
    .checkOvandsaln(ov, salnmask)
    ## get uniquely aligned-reads
    maskuniqaln <- .getMaskUniqueAln(alnreadids)
    if (ttpar@singleEnd == TRUE) {
        # unique + multi-mapping reads (only once) are considered
        avgreadlen <- mean(avgreadlen[!duplicated(alnreadids)])
    } else {
        # only unique reads are considered
        avgreadlen <- mean(avgreadlen[(d <= 500) & maskuniqaln])
    }
    ## fetch all different read identifiers from the overlapping alignments
    readids <- unique(alnreadids[queryHits(ov)])
    ## fetch all different transcripts from the overlapping alignments
    tx_idx <- sort(unique(subjectHits(ov)))
    
    cntvec <- .ttQuantExpress(ov, alnreadids, readids, tx_idx, ttpar, iste,
                                maskuniqaln, avgreadlen)
    ## Original TEtranscripts coerces fractional counts to integer
    setNames(as.integer(cntvec), names(cntvec))
}

#' @importFrom GenomicAlignments extractAlignmentRangesOnReference
#' @importFrom GenomicAlignments cigar start qwidth first second
.getAveLen <- function(ttpar, alnreads) {
    if (is(alnreads, "GAlignments")) {
        avgreadlenaln <- qwidth(alnreads)
    } else if (is(alnreads, "GAlignmentPairs")) {
        d <- abs(start(first(alnreads)) - start(second(alnreads)))
        d2 <- qwidth(second(alnreads))
        avgreadlenaln <- d + d2
    } else if (is(alnreads, "GAlignmentsList")) {
        d <- max(abs(diff(start(alnreads))))
        d[d<0] <- 0
        l <- lengths(alnreads)
        d2 <- unlist(qwidth(alnreads))[cumsum(l)]
        avgreadlenaln <- d + d2
    } else {
        stop(sprintf(".getAveLen: wrong class %s\n", class(alnreads)))
    }
    avgreadlenaln
}

.ttQuantExpress <- function(ov, alnreadids, readids, tx_idx, ttpar, iste,
                            maskuniqaln, avgreadlen) {
    
    ## build a matrix representation of the overlapping alignments
    ovalnmat <- .buildOvAlignmentsMatrix(ov, alnreadids, readids, tx_idx)
    
    mt <- match(readids, alnreadids)
    multigcnt <- rep(0L, length(ttpar@features))
    istex <- as.vector(iste[tx_idx])
    
    if (!all(iste)) {
        ## Correcting for preference of unique/multi-mapping reads to genes or TEs
        ovalnmat <- .correctPreference(ovalnmat, maskuniqaln, mt, istex)
    }
    
    # Getting counts from unique reads
    uniqcnt <- .countUniqueRead(ttpar, ovalnmat, maskuniqaln, mt, tx_idx, istex)
    
    if (!all(iste) & any(!maskuniqaln)) {
        ## Getting gene counts where a multimapping read maps to > 1 gene
        multigcnt <- .countMultiReadsGenes(ttpar, ovalnmat, maskuniqaln, mt,
                                        iste, istex, tx_idx, readids, 
                                        alnreadids, ov, uniqcnt)
    }
    
    # Adjusting 'ovalnmat' when alignment from multimapping read maps > 1 TE
    ovalnmat <- .adjustalnmultiovTE(iste, ov, alnreadids, readids, tx_idx, ovalnmat)
    
    cntvec <- rep(0L, length(ttpar@features))
    cntvec <- .ttEMstep(maskuniqaln, mt, ovalnmat, istex, tx_idx, readids, ttpar,
                        avgreadlen, cntvec)
    ## Summarize counts of unique and multimapping reads
    cntvec <- .summarizeCounts(iste, cntvec, uniqcnt, multigcnt, ttpar)
    cntvec
}



#' @importFrom Matrix Matrix sparseMatrix t which
#' @importFrom sparseMatrixStats colSums2 rowSums2
#' @importFrom SQUAREM squarem
#' @importFrom IRanges ranges
.ttEMstep <- function(maskuniqaln, mt, ovalnmat, istex, tx_idx, readids, ttpar,
                        avgreadlen, cntvec) {
if (sum(!maskuniqaln[mt]) > 0) { ## multi-mapping reads
    ## TEtranscripts doesn't use uniquely-aligned reads to inform the
    ## procedure of distributing multiple-mapping reads, as explained in
    ## Jin et al. (2015) pg. 3594, "to reduce potential bias to certain TEs."
    ## Hence, once counted, we discard unique alignments
    ovalnmat <- ovalnmat[!maskuniqaln[mt], istex]
    yesov <- rowSums2(ovalnmat)>0
    ovalnmat <- ovalnmat[yesov,]
    readids <- readids[!maskuniqaln[mt]][yesov]
    ## the Qmat matrix stores row-wise the probability that read i maps to
    ## a transcript j, assume uniform probabilities by now
    Qmat <- Matrix(0, nrow=length(readids), ncol=length(tx_idx[istex]),
                    dimnames=list(readids, NULL))
    Qmat[which(ovalnmat>0, arr.ind=TRUE)] <- 1
    # Qmat <- Qmat / rowSums2(ovalnmat)
    
    ## Pi, corresponding to rho in Equations (1), (2) and (3) in Jin et al.
    ## (2015) stores probabilities of expression for each transcript, corrected
    ## for its effective length as defined in Eq. (1) of Jin et al. (2015)
    # Pi <- colSums2(Qmat)
    # Pi <- colSums2(ovalnmat / rowSums2(ovalnmat))
    Pi <- colSums2(ovalnmat) / sum(ovalnmat)
    
    if (is(ttpar@features,"GRangesList")) {
        elen <- as.numeric(sum(width(ttpar@features[tx_idx][istex]))) - avgreadlen+1
    } else {
        elen <- width(ttpar@features[tx_idx][istex]) - avgreadlen + 1
    }
    Pi <- .correctForTxEffectiveLength(Pi, elen)
    
    ## use the SQUAREM algorithm to achieve faster EM convergence
    emres <- squarem(par=Pi, Q=Qmat, elen=elen,
                    fixptfn=.ttFixedPointFun,
                    control=list(tol=ttpar@tolerance, maxiter=ttpar@maxIter))
    Pi <- emres$par
    Pi[Pi < 0] <- 0 ## Pi estimates are sometimes negatively close to zero
    Pi <- Pi / sum(Pi)
    ## use the estimated transcript expression probabilities
    ## to finally distribute ambiguously mapping reads
    # probmassbyread <- as.vector(ovalnmat %*% Pi) 
    probmassbyread <- as.vector(Qmat %*% Pi)
    
    # Version 1
    # cntvecovtx <- rowSums(t(ovalnmat / probmassbyread) * Pi, na.rm=TRUE)
    # Version 2
    # wh <- which(ovalnmat, arr.ind=TRUE)
    cntvecovtx <- rep(0, ncol(ovalnmat)) #cntvecovtx <- rep(0, length(tx_idx))
    # x <- tapply(Pi[wh[, "col"]] / probmassbyread[wh[, "row"]], wh[, "col"],
    #             FUN=sum, na.rm=TRUE)
    wh <- which(Qmat==1, arr.ind=TRUE)
    x <- tapply(Pi[wh[, "col"]] / probmassbyread[wh[, "row"]], wh[, "col"],
                FUN=sum, na.rm=TRUE)
    
    cntvecovtx[as.integer(names(x))] <- x
    cntvec[tx_idx][istex] <- cntvecovtx
}
cntvec
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

#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors queryHits subjectHits
.buildOvAlignmentsMatrix <- function(ov, arids, rids, fidx) {
    mt1 <- match(arids[queryHits(ov)], rids)
    mt2 <- match(subjectHits(ov), fidx)
    positions <- cbind(mt1, mt2)
    oamat <- sparseMatrix(positions[,1], positions[,2], x=TRUE, 
                        dims = c(length(rids), length(fidx)))
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
    x[elen <= 0] <- 0
    if (sum(x) > 0) {
        x <- x / sum(x)
    }
    x
}

## private function .ttEstep()
## E-step of the EM algorithm of TEtranscripts
.ttEstep <- function(Q, Pi) {
    ## this first part is a sparse optimization of doing
    ## X <- t(t(Q) * Pi)
    X <- Q
    j <- rep(seq_len(ncol(X)), diff(X@p))
    X@x <- X@x * Pi[j]
    ## End of first part
    X <- X[rowSums2(X)>0,, drop=FALSE]
    X <- X / rowSums2(X)
    X
}

## private function .ttEstep()
## M-step of the EM algorithm of TEtranscripts
.ttMstep <- function(X) {
    Pi <- colSums2(X) / sum(X@x)
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


## private function .correctPreference()
## Corrects ovalnmat for preference of unique/multi-mapping reads to
## genes/TEs, respectively
.correctPreference <- function(ovalnmat, maskuniqaln, mt, istex) {
    indx <- (rowSums2(ovalnmat[,istex]) > 0) & (rowSums2(ovalnmat[,!istex]) > 0)
    
    ## Assigning unique reads mapping to both a TE and a gene as gene counts
    # Which unique reads overlap to both genes and TEs?
    idxu <- indx & maskuniqaln[mt]
    if (any(idxu)) {
        # ovalnmat[idxu,istex] <- FALSE
        whu <- which(ovalnmat[idxu,istex], arr.ind = TRUE)
        if (is(whu, "integer")) {
            whu <- as.matrix(data.frame(row = 1, col = whu))
        }
        whudf <- cbind(which(idxu)[whu[,"row"]], which(istex)[whu[,"col"]])
        ovalnmat[whudf] <- FALSE
    }
    
    ## Removing overlaps of multi-mapping reads to genes if at least one
    ## alignment of the read overlaps a TE
    idxm <- indx & !maskuniqaln[mt]
    if (any(idxm)) {
        # ovalnmat[idxm,!istex] <- FALSE
        whm <- which(ovalnmat[idxm,!istex], arr.ind = TRUE)
        if (is(whm, "integer")) {
            whm <- as.matrix(data.frame(row = 1, col = whm))
        }
        whmdf <- cbind(which(idxm)[whm[,"row"]], which(!istex)[whm[,"col"]])
        ovalnmat[whmdf] <- FALSE
    }
    
    ovalnmat
}


## private function .countMultiReadsGenes()
## Counts multi-mapping reads mapping to multiple genes by counting fraction
## counts
.countMultiReadsGenes <- function(ttpar, ovalnmat, maskuniqaln, mt, iste,
                                    istex, tx_idx, readids, alnreadids, ov,
                                    uniqcnt) {
    ovalnmat_multig <- ovalnmat[!maskuniqaln[mt], !istex, drop=FALSE]
    yesg <- rowSums2(ovalnmat_multig)>0
    
    ## Computing counts for reads with multiple alignments mapping to different
    ## reads and also for multimapping reads with only 1 alignment mapping to
    ## a gene
    ovalnmat_multig <- ovalnmat_multig[yesg,,drop=FALSE]
    
    # Getting the num of different alignments mapping to a gene for each read
    alnreadids_multig <-alnreadids[unique(queryHits(
                                            ov[!iste[subjectHits(ov)]]))]
    nalnperread <- table(alnreadids_multig) # getting only overlaps from genes
    readids_multig <- readids[!maskuniqaln[mt]][yesg]
    mt_multig <- match(readids_multig, names(nalnperread))
    
    # Counts provided by unique reads for genes to which multimapping 
    # reads map to
    matmultiguniqc <- t(t(ovalnmat_multig)*uniqcnt[tx_idx][!istex])
    rsum <- rowSums2(matmultiguniqc)
    rsum0 <- which(rsum == 0)
    
    # Reads for which the genes to which the read aligns have > 0 counts
    # matmultiguniqc[!rsum0,, drop=FALSE] <- matmultiguniqc[!rsum0,,drop=FALSE]/
    #     rsum[!rsum0]
    j <- matmultiguniqc@i + 1
    matmultiguniqc@x[!(j %in% rsum0)] <- matmultiguniqc@x[!(j %in% rsum0)] / 
        rsum[j[!(j %in% rsum0)]]
    
    
    # Reads for which the genes to which the read aligns have 0 counts
    # matmultiguniqc[rsum0,, drop=FALSE] <- (ovalnmat_multig[rsum0,, drop=FALSE] /
    #     rowSums2(ovalnmat_multig[rsum0,, drop=FALSE]))
    
    matmultiguniqc@x[(j %in% rsum0)] <- ovalnmat_multig@x[(j %in% rsum0)] / 
        rowSums2(ovalnmat_multig)[j[j %in% rsum0]]
    
    # Adjusting for number of alignments
    matmultiguniqc <- matmultiguniqc/as.numeric(nalnperread[mt_multig])
    
    multigcnt <- rep(0L, length(ttpar@features))
    multigcnt[tx_idx][!istex] <- colSums2(matmultiguniqc)
    multigcnt
}

## private function .countUniqueRead()
## Counts unique reads mapping to TEs and genes (if present)
.countUniqueRead <- function(ttpar, ovalnmat, maskuniqaln, mt, tx_idx, istex) {
    uniqcnt <- rep(0L, length(ttpar@features))
    ovalnmatuniq_g <- ovalnmat[maskuniqaln[mt], !istex, drop=FALSE]
    ovalnmatuniq_te <- ovalnmat[maskuniqaln[mt], istex, drop=FALSE]
    ovmultig <- rowSums2(ovalnmatuniq_g) > 1
    ovmultite <- rowSums2(ovalnmatuniq_te) > 1
    
    # Addressing gene counts: count divided by the number of different genes
    if (any(ovmultig)) {
        # Counting reads overlapping only 1 element
        uniqcnt[tx_idx][!istex] <- colSums2(ovalnmatuniq_g[!ovmultig,,drop=FALSE])
        # Counting reads overlapping more than 1 element
        mg <- ovalnmatuniq_g[ovmultig,,drop=FALSE] /
            rowSums2(ovalnmatuniq_g[ovmultig,,drop=FALSE])
        uniqcnt[tx_idx][!istex] <- uniqcnt[tx_idx][!istex] + colSums2(mg)
    
    } else {
        uniqcnt[tx_idx][!istex] <- colSums2(ovalnmatuniq_g)
    }
    
    # Addressing TE counts: count divided by the number of different TEs
    # proportionally to the expression level of each TE provided by unique counts
    if (any(ovmultite)) {
        # Counting reads overlapping only 1 element
        uniqcnt[tx_idx][istex] <- colSums2(ovalnmatuniq_te[!ovmultite,,drop=FALSE])
        # Counting reads overlapping more than 1 element
        mte <- t(t(
            ovalnmatuniq_te[ovmultite,,drop=FALSE])*uniqcnt[tx_idx][istex])
        nocounts <- which(rowSums2(mte) == 0)
        mte[nocounts,,drop=FALSE] <- ovalnmatuniq_te[
            ovmultite,,drop=FALSE][nocounts,]
        mte <- mte/rowSums2(mte)
        uniqcnt[tx_idx][istex] <- uniqcnt[tx_idx][istex] + colSums2(mte)
    
    } else {
        uniqcnt[tx_idx][istex] <- colSums2(ovalnmatuniq_te)
    }
    
    uniqcnt
}


## private function .summarizeCounts()
## Counts unique reads mapping to TEs and genes (if present)
.summarizeCounts <- function(iste, cntvec, uniqcnt, multigcnt, ttpar) {
    
    ## add multi-mapping and unique-mapping counts
    if (!all(iste)) {
        cntvec <- cntvec + uniqcnt + multigcnt
    } else {
        cntvec <- cntvec + uniqcnt
    }
    names(cntvec) <- names(ttpar@features)
    
    cntvec_t <- cntvec[iste]
    ## aggregate TE quantifications if necessary
    if (length(ttpar@aggregateby) > 0) {
        f <- .factoraggregateby(ttpar@features[iste], ttpar@aggregateby)
        stopifnot(length(f) == length(cntvec_t)) ## QC
        cntvec_t <- tapply(cntvec_t, f, sum, na.rm=TRUE)
    }
    
    ## aggregating exon counts to genes
    if (!all(iste)) {
        cntvec <- c(cntvec_t, cntvec[!iste])
    } else {
        cntvec <- cntvec_t
    }
    
    cntvec
}

## private function .adjustalnmultiovTE()
## Adjusts ovalnmat when a multimapping reads have alignments mapping to > 1
## TE, by dividing the count value by the number of different TEs the
## alignment maps to

#' @importFrom S4Vectors queryHits subjectHits
.adjustalnmultiovTE <- function(iste, ov, alnreadids, readids, tx_idx,
                                ovalnmat) {
    isteov <- iste[subjectHits(ov)]
    multialn <- !(.getMaskUniqueAln(queryHits(ov)))
    ovte_multialn <- ov[isteov & multialn]
    mt1_tem <- match(alnreadids[queryHits(ovte_multialn)], readids)
    mt2_tem <- match(subjectHits(ovte_multialn), tx_idx)
    pos_tem <- cbind(mt1_tem, mt2_tem)
    novperaln <- table(queryHits(ovte_multialn))
    mt_multig <- match(queryHits(ovte_multialn), names(novperaln))
    ovalnmat <- as(as(as(ovalnmat, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    ovalnmat[pos_tem] <- ovalnmat[pos_tem] / as.numeric(novperaln[mt_multig])
    ovalnmat
}
