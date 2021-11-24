#' Build a Telescope parameter object
#'
#' Build an object of the class \code{TelescopeParam}.
#'
#' @param bfl A \code{BamFile} or \code{BamFileList} object, or a character
#' string vector of BAM filenames.
#'
#' @param teFeatures A \code{GRanges} or \code{GRangesList} object. Elements
#' in this object should have names, which will be used as a grouping factor
#' for ranges forming a common locus (equivalent to "locus" column in
#' Telescope), unless other metadata column names are specified in the
#' \code{aggregateby} parameter.
#' 
#' @param aggregateby Character vector with column names from the annotation
#' to be used to aggregate quantifications. By default, this is an empty
#' vector, which means that the names of the input \code{GRanges} or
#' \code{GRangesList} object given in the \code{teFeatures} parameter are used
#' to aggregate quantifications.
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
#' When \code{fragments=FALSE}, the read-counting method only counts
#' ‘mated pairs’ from opposite strands, while when \code{fragments=TRUE}
#' (default), same-strand pairs, singletons, reads with unmapped pairs and
#' other fragments are also counted. \code{fragments=TRUE} is equivalent to
#' the original Telescope algorithm. For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}()}.
#' 
#' @param minOverlFract (Default 0.2) A numeric scalar. \code{minOverlFract}
#' is multiplied by the median read length and the resulting value is used to
#' specify the \code{minoverlap} argument from
#' \code{\link[IRanges:findOverlaps-methods]{findOverlaps}} from the
#' \pkg{IRanges} package. When no minimum overlap is required, set
#' \code{minOverlFract = 0}.
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
#' TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
#'                     package="atena"))
#' gene_annot <- readRDS(file = system.file("extdata", "Top50genes.rds",
#'                                          package="atena"))
#' tspar <- TelescopeParam(bfl=bamfiles, teFeatures=TE_annot,
#'                         geneFeatures = gene_annot,
#'                         singleEnd = TRUE, ignoreStrand=TRUE)
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
#' @rdname TelescopeParam-class
TelescopeParam <- function(bfl, teFeatures, aggregateby=character(0),
                            geneFeatures=NA,
                            singleEnd=TRUE, 
                            strandMode=1L, 
                            ignoreStrand=FALSE,
                            fragments=FALSE,
                            minOverlFract=0.2,
                            pi_prior=0L, 
                            theta_prior=0L, 
                            em_epsilon=1e-7,
                            maxIter=100L) {
    bfl <- .checkBamFileListArgs(bfl, singleEnd, fragments)
    
    features <- .processFeatures(teFeatures, deparse(substitute(teFeatures)),
                                geneFeatures,deparse(substitute(geneFeatures)),
                                aggregateby, aggregateexons = TRUE)
    
    new("TelescopeParam", bfl=bfl, features=features,
        aggregateby=aggregateby, singleEnd=singleEnd,ignoreStrand=ignoreStrand,
        strandMode=as.integer(strandMode), fragments=fragments,
        minOverlFract=minOverlFract, pi_prior=pi_prior,
        theta_prior=theta_prior, em_epsilon=em_epsilon,
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
            
            features <- .consolidateFeatures(x, rownames(cnt)[-nrow(cnt)])
            SummarizedExperiment(assays=list(counts=cnt),
                                    rowRanges=c(features),
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
    asvalues <- integer()
    mreadlen <- numeric(0)
    strand_arg <- "strandMode" %in% formalArgs(readfun)
    yieldSize(bf) <- yieldSize
    open(bf)
    while (length(alnreads <- do.call(readfun,
                                c(list(file = bf), list(param=param),
                                list(strandMode=tspar@strandMode)[strand_arg],
                                list(use.names=TRUE))))) {
        alnreadids <- c(alnreadids, names(alnreads))
        asvalues <- c(asvalues, mcols(alnreads)$AS)
        mreadlen <- median(width(ranges(alnreads)))
        thisov <- mode(alnreads, tspar@features,
                        minOverlFract=as.integer(tspar@minOverlFract*mreadlen),
                        ignoreStrand=tspar@ignoreStrand)
        ov <- .appendHits(ov, thisov)
    }
    # close(bf)
    on.exit(close(bf))
    maskuniqaln <- !(duplicated(alnreadids) | 
                        duplicated(alnreadids, fromLast = TRUE))
    ## fetch all different read identifiers from the overlapping alignments
    readids <- unique(alnreadids[queryHits(ov)])
    ## Adding "no_feature" overlaps to 'ov'
    ov <- .getNoFeatureOv(maskuniqaln, ov, alnreadids)
    mt <- match(readids, alnreadids)
    readids <- unique(alnreadids[queryHits(ov)]) # updating 'readids'
    cntvec <- .tsEMstep(tspar, alnreadids, readids, ov, asvalues,
                        iste, maskuniqaln, mt)
    setNames(as.integer(cntvec), names(cntvec))
}

## private function .getNoFeatureOv(), obtains the overlaps to "no_feature"
#' @importFrom S4Vectors Hits queryHits subjectHits nRnode nLnode from to
.getNoFeatureOv <- function(maskuniqaln, ov, alnreadids) {
    maskmultialn_ov <- !(maskuniqaln[queryHits(ov)])
    alnreadidx_all <- match(alnreadids, unique(alnreadids))
    ovid <- unique(alnreadidx_all[unique(queryHits(ov)[maskmultialn_ov])])
    alnall <- which(alnreadidx_all %in% ovid)
    nofeat <- setdiff(alnall, queryHits(ov)[maskmultialn_ov])
    from <- nofeat
    to <- rep(nRnode(ov) + 1, length(nofeat))
    nofeat_hits <- Hits(from = from, to = to, nLnode = nLnode(ov), 
                        nRnode = max(to), sort.by.query = TRUE)
    
    # Adding overlaps to "no_feature" to 'ov'
    hits1 <- ov
    hits2 <- nofeat_hits
    
    hits <- c(Hits(from=from(hits1), to=to(hits1),
                    nLnode=nLnode(hits1),
                    nRnode=nRnode(hits2), sort.by.query=FALSE),
                Hits(from=from(hits2), to=to(hits2),
                    nLnode=nLnode(hits2),
                    nRnode=nRnode(hits2), sort.by.query=FALSE))
    
    ovnofeat <- Hits(from = from(hits), to = to(hits), nLnode = nLnode(hits),
                    nRnode = nRnode(hits), sort.by.query=TRUE)
    return(ovnofeat)
}

#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom Matrix Matrix rowSums colSums t which
#' @importFrom SQUAREM squarem
.tsEMstep <- function(tspar, alnreadids, readids, ov, asvalues,
                      iste, maskuniqaln, mt) {
    ## initialize vector of counts derived from multi-mapping reads
    cntvec <- rep(0L, length(tspar@features) + 1)
    
    # Getting counts from reads overlapping only one feature (this excludes --> SOBRA
    # unique reads mapping to two or more overlapping features)
    # ovunique <- rowSums(ovalnmat) == 1
    # cntvec[tx_idx] <- colSums(ovalnmat[ovunique,])
    
    alnreadidx <- match(alnreadids, readids)
    rd_idx <- sort(unique(alnreadidx[queryHits(ov)]))
    
    ## fetch all different transcripts from the overlapping alignments
    tx_idx <- sort(unique(subjectHits(ov)))
    istex <- as.vector(iste[tx_idx])[-length(tx_idx)]  # removing "no_feature"
    
    asvalues <- (asvalues-min(asvalues)+1) / (max(asvalues)+1 - min(asvalues))
    QmatTS <- .buildOvValuesMatrix(tspar, ov, asvalues, alnreadidx,
                                    rd_idx, tx_idx)
    
    if (!all(iste)) {
        ## Correcting for preference of unique/multimapping reads to genes/TEs
        QmatTS_nofeat <- QmatTS[,ncol(QmatTS), drop = FALSE]
        QmatTS <- .correctPreferenceTS(QmatTS[,-ncol(QmatTS)], maskuniqaln,
                                       mt, istex)
        stopifnot(!any(rowSums(QmatTS[,istex]) > 0 &
                           rowSums(QmatTS[,!istex]) > 0))
        QmatTS <- cbind(QmatTS, QmatTS_nofeat)
    }
    
    # --- EM-step --- 
    QmatTS <- QmatTS / rowSums(QmatTS)
    
    # Getting 'y' indicator of unique and multi-mapping status from QmatTS 
    maskmulti <- ifelse(rowSums(QmatTS > 0) == 1, 0, 1)
    
    # Telescope (Bendall et al.(2019)) defines the initial π estimate uniformly
    PiTS <- rep(1 / length(tx_idx), length(tx_idx))
    
    # Telescope defines an additional reassignment parameter (θ) as uniform
    Theta <- rep(1 / length(tx_idx), length(tx_idx))
    
    ## The SQUAREM algorithm to run the EM procedure
    a <- tspar@pi_prior # 0
    b <- tspar@theta_prior # 0
    Thetaenv <- new.env()
    assign("Theta", Theta, envir=Thetaenv)
    tsres <- squarem(par=PiTS, Thetaenv=Thetaenv, Q=QmatTS,maskmulti=maskmulti,
                    a=a, b=b, fixptfn=.tsFixedPointFun,
                    control=list(tol=tspar@em_epsilon, maxiter=tspar@maxIter))
    PiTS <- tsres$par
    PiTS[PiTS < 0] <- 0 ## Pi estimates are sometimes negatively close to zero
    # --- end EM-step ---
    
    # Implementation for reassign_mode=exclude
    Theta <- get("Theta", envir=Thetaenv)
    X <- .tsEstep(QmatTS, Theta, maskmulti, PiTS)
    maxbyrow <- rowMaxs(X)
    # Xind <- X == maxbyrow
    # Quicker version of the previous line
    Xind <- (X / maxbyrow) == 1
    nmaxbyrow <- rowSums(Xind)
    
    # Do not differenciate between TE and genes with istex because in Telescope
    # genes are also included in the EMstep.
    cntvec[tx_idx] <- colSums(Xind[nmaxbyrow == 1, ])
    
    names(cntvec) <- c(names(tspar@features), "no_feature")
    nofeat <- cntvec["no_feature"]
    cntvec <- .tssummarizeCounts(cntvec[-length(cntvec)], iste, tspar)
    cntvec <- c(cntvec, nofeat)
    cntvec
}

## private function .correctPreferenceTS()
## Corrects QmatTS for preference of unique/multi-mapping reads to
## genes/TEs, respectively, in Telescope
.correctPreferenceTS <- function(QmatTS, maskuniqaln, mt, istex) {
    indx <- (rowSums(QmatTS[,istex]) > 0) & (rowSums(QmatTS[,!istex]) > 0)
    
    ## Assigning unique reads mapping to both a TE and a gene as gene counts
    # Which unique reads overlap to both genes and TEs?
    idxu <- indx & maskuniqaln[mt]
    if (any(idxu)) {
        # QmatTS[idxu,istex] <- FALSE
        whu <- which(QmatTS[idxu,istex] > 0, arr.ind = TRUE)
        whudf <- cbind(which(idxu)[whu[,"row"]], which(istex)[whu[,"col"]])
        QmatTS[whudf] <- FALSE
    }
    
    ## Removing overlaps of multi-mapping reads to genes if at least one
    ## alignment of the read overlaps a TE
    idxm <- indx & !maskuniqaln[mt]
    if (any(idxm)) {
        # QmatTS[idxm,!istex] <- FALSE
        whm <- which(QmatTS[idxm,!istex] > 0, arr.ind = TRUE)
        whmdf <- cbind(which(idxm)[whm[,"row"]], which(!istex)[whm[,"col"]])
        QmatTS[whmdf] <- FALSE
    }
    
    QmatTS
}

.tssummarizeCounts <- function(cntvec, iste, tspar) {
    cntvec_t <- cntvec[iste]
    ## aggregate TE quantifications if necessary
    if (length(tspar@aggregateby) > 0) {
        f <- .factoraggregateby(tspar@features[iste], tspar@aggregateby)
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




## private function .tsEstep()
## E-step of the EM algorithm of Telescope
.tsEstep <- function(Q, Theta, maskmulti, Pi) {
    X <- t(t(Q) * Pi)
    # X[maskmulti, ] <- t(t(X[maskmulti, ]) * Theta)
    # quicker computation of previous line
    wh <- which(X*maskmulti > 0, arr.ind = TRUE)
    X[wh] <- t(t(X[wh]) * Theta[wh[, "col"]])
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
    Theta <- ((colSums(X[maskmulti, , drop=FALSE]) + b) / 
        (sum(maskmulti) + b*ncol(X)))
}

## private function .tsFixedPointFun()
## fixed point function of the EM algorithm of Telescope
.tsFixedPointFun <- function(Pi, Thetaenv, Q, maskmulti, a, b) {
    Theta <- get("Theta", envir=Thetaenv)
    X <- .tsEstep(Q, Theta, maskmulti, Pi)
    Pi2 <- .tsMstepPi(X, a)
    Theta2 <- .tsMstepTheta (X, maskmulti, b)
    assign("Theta", Theta2, envir=Thetaenv)
    
    Pi2
}

