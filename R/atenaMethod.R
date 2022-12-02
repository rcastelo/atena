#' Build an atena parameter object
#'
#' Build an object of the class \code{atenaParam}.
#'
#' @param bfl A \code{BamFile} or \code{BamFileList} object, or a character
#' string vector of BAM filenames.
#'
#' @param teFeatures A \code{GRanges} or \code{GRangesList} object. Elements
#' in this object should have names, which are used as a grouping factor for
#' genomic ranges forming a common locus. This grouping is performed previous
#' to TE expression quantification, unlike the aggregation of quantifications 
#' performed when the \code{aggregateby} parameter is specified, which is 
#' performed after individual TE instances are quantified.
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
#' Ambiguous alignments (alignments overlapping > 1 feature) are not counted.
#'
#' @param geneFeatures A \code{GRanges} or \code{GRangesList} object with the
#' gene annotated features to be quantified. Unique reads are first tallied
#' with respect to these gene features whereas multi-mapping reads are
#' preferentially assigned to TEs. Elements should have names indicating the
#' gene name/id. In case that \code{geneFeatures} contains a metadata column
#' named \code{type}, only the elements with \code{type} = \code{exon} are 
#' considered for the analysis. Then, exon counts are summarized to the gene
#' level.
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
#' @param pi_prior (Default 0) A positive numeric object indicating the prior
#' on pi. The same prior can be specified for all features setting
#' \code{pi_prior} as a scalar, or each feature can have a specific prior by
#' setting \code{pi_prior} as a vector with \code{names()} corresponding to
#' all feature names. Setting a pi prior is equivalent to adding n unique
#' reads. 
#'
#' @param theta_prior (Default 0) A positive numeric object indicating the
#' prior on Q. The same prior can be specified for all features setting
#' \code{theta_prior} as a scalar, or each feature can have a specific prior by
#' setting \code{theta_prior} as a vector with \code{names()} corresponding to
#' all feature names. Equivalent to adding n non-unique reads.
#'
#' @param em_epsilon (Default 1e-7) A numeric scalar indicating the EM
#' Algorithm Epsilon cutoff.
#'
#' @param maxIter A positive integer scalar storing the maximum number of
#' iterations of the EM SQUAREM algorithm (Du and Varadhan, 2020). Default
#' is 100 and this value is passed to the \code{maxiter} parameter of the
#' \code{\link[SQUAREM]{squarem}()} function.
#' 
#' @param reassign_mode (Default 'exclude') Character vector indicating
#' reassignment mode after EM step. 
#' Available methods are 'exclude' (reads with more than one best
#' assignment are excluded from the final counts), 'choose' (when reads have
#' more than one best assignment, one of them is randomly chosen), 'average'
#' (the read count is divided evenly among the best assignments) and 'conf'
#' (only assignments that exceed a certain threshold -defined by 
#' \code{conf_prob} parameter- are accepted, then the read count is
#' proportionally divided among the assignments above \code{conf_prob}).
#' 
#' @param conf_prob (Default 0.9) Minimum probability for high confidence
#' assignment.
#'
#'
#' @details
#' This is the constructor function for objects of the class
#' \code{atenaParam-class}. This type of object is the input to the
#' function \code{\link{qtex}()} for quantifying expression of transposable
#' elements, which will call the atena method with this type of object. The
#' atena method uses a multiple '__no_feature' approach in which as many
#' '__no_feature' features as different overlapping patterns of multimapping
#' reads in the overlapping matrix are used to represent alignments mapping
#' outside annotations.
#'
#' @return A \linkS4class{atenaParam} object.
#'
#' @examples
#' bamfiles <- list.files(system.file("extdata", package="atena"),
#'                        pattern="*.bam", full.names=TRUE)
#' TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
#'                     package="atena"))
#' gene_annot <- readRDS(file = system.file("extdata", "Top50genes.rds",
#'                                          package="atena"))
#' atpar <- atenaParam(bfl=bamfiles, teFeatures=TE_annot,
#'                         geneFeatures = gene_annot,
#'                         singleEnd = TRUE, ignoreStrand=TRUE)
#' atpar
#'
#'
#' @importFrom methods is new
#' @importFrom Rsamtools BamFileList
#' @importFrom S4Vectors mcols
#' @export
#' @rdname atenaParam-class
atenaParam <- function(bfl, teFeatures, aggregateby=character(0),
                            ovMode="ovUnion",
                            geneFeatures=NA,
                            singleEnd=TRUE,
                            strandMode=1L,
                            ignoreStrand=FALSE,
                            fragments=FALSE,
                            pi_prior=0L,
                            theta_prior=0L,
                            em_epsilon=1e-7,
                            maxIter=100L,
                            reassign_mode="exclude",
                            conf_prob=0.9) {
    bfl <- .checkBamFileListArgs(bfl, singleEnd, fragments)
    
    if (!reassign_mode %in% c("exclude","choose","average","conf"))
      stop("'reassign_mode' should be one of 'exclude', 'choose', 'average' or 'conf'")
    
    if (!ovMode %in% c("ovUnion","ovIntersectionStrict"))
      stop("'ovMode' should be one of 'ovUnion', 'ovIntersectionStrict'")
    
    features <- .processFeatures(teFeatures, deparse(substitute(teFeatures)),
                                geneFeatures,deparse(substitute(geneFeatures)),
                                aggregateby, aggregateexons = TRUE)
    .checkPriors(names(features), names(pi_prior), names(theta_prior))
    
    new("atenaParam", bfl=bfl, features=features,
        aggregateby=aggregateby, ovMode=ovMode,
        singleEnd=singleEnd,ignoreStrand=ignoreStrand,
        strandMode=as.integer(strandMode), fragments=fragments,
        pi_prior=pi_prior,
        theta_prior=theta_prior, em_epsilon=em_epsilon,
        maxIter=as.integer(maxIter), reassign_mode=reassign_mode,
        conf_prob=conf_prob)
}

#' @param object A \linkS4class{atenaParam} object.
#'
#' @importFrom GenomeInfoDb seqlevels
#' @export
#' @aliases show,atenaParam-method
#' @rdname atenaParam-class
setMethod("show", "atenaParam",
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
#' @aliases qtex,atenaParam-method
#' @rdname qtex
setMethod("qtex", "atenaParam",
        function(x, phenodata=NULL, mode=ovUnion, yieldSize=1e6L,
                    BPPARAM=SerialParam(progressbar=TRUE)) {
            .checkPhenodata(phenodata, length(x@bfl))
            
            cnt <- bplapply(x@bfl, .qtex_atena, atpar=x, mode=x@ovMode,
                            yieldSize=yieldSize, BPPARAM=BPPARAM)
            cnt <- .processMultiNofeature(cnt, x)
            cnt <- do.call("cbind", cnt)
            colData <- .createColumnData(cnt, phenodata)
            colnames(cnt) <- rownames(colData)
            whnofeat <- grep(x = rownames(cnt), pattern = "no_feature")
            features <- .consolidateFeatures(x, rownames(cnt)[-whnofeat],
                                             whnofeat)
            # features <- .consolidateFeatures(x, rownames(cnt)[-nrow(cnt)])
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
#' @importFrom Matrix Matrix t which
#' @importFrom sparseMatrixStats rowSums2 colSums2
#' @importFrom SQUAREM squarem
.qtex_atena <- function(bf, atpar, mode, yieldSize=1e6L) {
    mode <- match.fun(mode)
    readfun <- .getReadFunction(atpar@singleEnd, atpar@fragments)
    .checkreassignModes(atpar)
    sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                            isDuplicate=FALSE,
                            isNotPassingQualityControls=FALSE)
    param <- ScanBamParam(flag=sbflags, what="flag", tag="AS")
    iste <- as.vector(attributes(atpar@features)$isTE[,1])
    if (any(duplicated(names(atpar@features[iste])))) {
        stop(".qtex_atena: transposable element annotations do not contain unique names for each element")
    }
    ov <- Hits(nLnode=0, nRnode=length(atpar@features), sort.by.query=TRUE)
    alnreadids <- character(0)
    asvalues <- integer()
    alen <- numeric(0)
    mreadlen <- numeric(0)
    salnmask <- logical(0)
    strand_arg <- "strandMode" %in% formalArgs(readfun)
    yieldSize(bf) <- yieldSize
    open(bf)
    while (length(alnreads <- do.call(readfun,
                                c(list(file = bf), list(param=param),
                                list(strandMode=atpar@strandMode)[strand_arg],
                                list(use.names=TRUE))))) {
        alnreadids <- c(alnreadids, names(alnreads))
        asvalues <- c(asvalues, .getAlignmentASScoreTS(alnreads, tag = "AS"))
        readlen <- .getAlignmentLength(alnreads)
        alen <- c(alen, readlen)
        salnmask <- c(salnmask, any(.secondaryAlignmentMask(alnreads)))
        thisov <- mode(alnreads, atpar@features,
                       ignoreStrand=atpar@ignoreStrand, inter.feature=TRUE)
        ov <- .appendHits(ov, thisov)
    }
    # close(bf)
    on.exit(close(bf))
    .checkOvandsaln(ov, salnmask)
    maskuniqaln <- .getMaskUniqueAln(alnreadids)
    ## fetch all different read identifiers from the overlapping alignments
    readids <- unique(alnreadids[queryHits(ov)])
    ## Adding "no_feature" overlaps to 'ov'
    ov <- .getNoFeatureOv(maskuniqaln, ov, alnreadids)
    mt <- match(readids, alnreadids)
    readids <- unique(alnreadids[queryHits(ov)]) # updating 'readids'
    cntvec <- .atEMstep(atpar, alnreadids, readids, ov, asvalues,
                        iste, maskuniqaln, mt, alen)
    setNames(as.integer(cntvec), names(cntvec))
}


#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom Matrix Matrix t which
#' @importFrom SQUAREM squarem
#' @importClassesFrom Matrix lgCMatrix
.atEMstep <- function(atpar, alnreadids, readids, ov, asvalues,
                      iste, maskuniqaln, mt, alen) {
    ## initialize vector of counts derived from multi-mapping reads
    cntvec <- rep(0L, length(atpar@features) + 1)
    
    alnreadidx <- match(alnreadids, readids)
    rd_idx <- sort(unique(alnreadidx[queryHits(ov)]))
    
    ## fetch all different transcripts from the overlapping alignments
    tx_idx <- sort(unique(subjectHits(ov)))
    if (all(maskuniqaln)) {
        istex <- as.vector(iste[tx_idx])
    } else {
        istex <- as.vector(iste[tx_idx])[-length(tx_idx)] # removing no_feature
    }
    # asvalues <- (asvalues-min(asvalues)+1) / (max(asvalues)+1 - min(asvalues))
    asvalues <- .rescaleASat(asvalues, alen = alen)
    QmatAT <- .buildOvValuesMatrix(atpar, ov, asvalues, alnreadidx,
                                    rd_idx, tx_idx)
    
    if (!all(iste)) { ## Correcting for preference of reads to genes/TEs
        if (all(maskuniqaln)) {
          QmatAT <- .correctPreferenceTS(QmatAT, maskuniqaln, mt, istex)
          stopifnot(!any(rowSums2(QmatAT[,istex]) > 0 &
                           rowSums2(QmatAT[,!istex]) > 0))
        } else {
            QmatAT_nofeat <- QmatAT[,ncol(QmatAT), drop = FALSE]
            QmatAT <- .correctPreferenceTS(QmatAT[,-ncol(QmatAT)], maskuniqaln,
                                           mt, istex)
            stopifnot(!any(rowSums2(QmatAT[,istex]) > 0 &
                               rowSums2(QmatAT[,!istex]) > 0))
            QmatAT <- cbind(QmatAT, QmatAT_nofeat)
        }
    }
    
    nfeatures <- ncol(QmatAT)
    QmatAT <- .multiNofeature(QmatAT)
    nnofeat <- ncol(QmatAT) - nfeatures + 1
    cntvec <- c(cntvec, rep(0L, nnofeat - 1))
    tx_idx <- c(tx_idx, (max(tx_idx)+1):(max(tx_idx)+nnofeat -1))
    
    # --- EM-step --- 
    # QmatAT <- QmatAT / rowSums2(QmatAT)
    QmatAT@x <- QmatAT@x / rowSums2(QmatAT)[QmatAT@i +1]
    
    # Getting 'y' indicator of unique and multi-mapping status from QmatAT 
    # maskmulti <- ifelse(rowSums2(QmatAT > 0) == 1, 0, 1)
    maskmulti <- rowSums2(QmatAT > 0) > 1
    
    # π estimate proportional to relative abundances of features
    Pi <- colSums2(QmatAT / rowSums2(QmatAT))
    
    # θ proportional to relative abundances of features given by multimapping reads
    Theta <- colSums2(QmatAT[maskmulti,] / rowSums2(QmatAT[maskmulti,]))
    
    ## The SQUAREM algorithm to run the EM procedure
    a <- .processPriors(atpar@pi_prior, atpar, tx_idx, nnofeat) 
    b <- .processPriors(atpar@theta_prior, atpar, tx_idx, nnofeat) 

    Thetaenv <- new.env()
    assign("Theta", Theta, envir=Thetaenv)
    atres <- squarem(par=Pi, Thetaenv=Thetaenv, Q=QmatAT,maskmulti=maskmulti,
                    a=a, b=b, fixptfn=.tsFixedPointFun,
                    control=list(tol=atpar@em_epsilon, maxiter=atpar@maxIter))
    Pi <- atres$par
    Pi[Pi < 0] <- 0 ## Pi estimates are sometimes negatively close to zero
    # --- end EM-step ---
    
    Theta <- get("Theta", envir=Thetaenv)
    X <- .tsEstep(QmatAT, Theta, maskmulti, Pi)
    cntvec <- .reassign(X, atpar@reassign_mode, atpar@conf_prob, cntvec, tx_idx)
    
    nofeat_names <- paste("no_feature", 1:nnofeat, sep = "")
    names(cntvec) <- c(names(atpar@features), nofeat_names)
    nofeat <- cntvec[nofeat_names]
    cntvec <- .tssummarizeCounts(cntvec[1:length(atpar@features)], iste, atpar)
    cntvec <- c(cntvec, nofeat)
    cntvec
}


## private function .multiNofeature()
## Implements multiple '__no_feature' method
#' @importFrom stats aggregate
#' @importFrom Matrix Matrix
.multiNofeature <- function(QmatAT) { 
  ## Implementation of multiple no_feature method
  ovalnmat <- QmatAT > 0
  comb <-  aggregate(summary(ovalnmat)$j, by = list(ovalnmat@i + 1), 
                     FUN = function(x) paste0(x, collapse = "_"))
  combid <- match(comb$x, unique(comb$x)) # col id
  QmatAT[,ncol(QmatAT)] # values no_feature
  1:nrow(QmatAT) # row id
  ovmatnofeat <- Matrix(do.call("numeric", list(1)),
                        nrow=nrow(QmatAT), ncol=max(combid))
  ovmatnofeat[cbind(1:nrow(QmatAT), combid)] <- QmatAT[,ncol(QmatAT)]
  
  QmatAT_nofeatmulti <- cbind(QmatAT[,-ncol(QmatAT)], ovmatnofeat)
  return(QmatAT_nofeatmulti)
}

## private function .processMultiNofeature()
## Addresses different number of '__no_feature' features from the multiple
## '__no_feature' method
.processMultiNofeature <- function(cnt, x) {
  # adapted from cbindX() from gdata package
  if (length(cnt) > 1) {
    nrowsam <- lengths(cnt)
    maxi <- which.max(nrowsam)
    test <- nrowsam < nrowsam[maxi]
    for(i in 1:length(nrowsam)) {
      if(test[i]) {
        add <- rep(0, nrowsam[maxi] - nrowsam[i])
        cnt[[i]] <- c(cnt[[i]], add)
        names(cnt[[i]]) <- names(cnt[[maxi]])
      }
    }
  }
  return(cnt)
}



#' @importFrom matrixStats logSumExp
.rescaleASat <- function(asvalues, alen) {
  # asrescale <- rescale(asvalues, c(1, max(asvalues) - min(asvalues) + 1))
  # asrescale <- asrescale + alen
  # asrescale <- asrescale/max(asrescale)
  # asrescale <- expm1(asrescale*100)
  asrescale <- exp(asvalues - logSumExp(asvalues))
  asrescale
}


## private function .checkPriors()
## Checks that if length(prior) > 1, all features after .processFeatures()
## have a corresponding prior with the same name
.checkPriors <- function(namesfeatures, namespriora, namespriorb) {
  if (length(namespriora) > 1 & (!(all(namesfeatures %in% namespriora))))
    stop("after processing annotations (.processFeatures()), not all features have a pi_prior matching its name. Check that 'names()' of 'teFeatures' or 'geneFeatures' match those of 'pi_prior'")
  
  if (length(namespriorb) > 1 & (!(all(namesfeatures %in% namespriorb))))
    stop("after processing annotations (.processFeatures()), not all features have a theta_prior matching its name. Check that 'names()' of 'teFeatures' or 'geneFeatures' match those of 'theta_prior'")
}


## private function .processPriors()
## Assigns priors to "no_feature" and checks the length of the vector
## with priors
.processPriors <- function(prior, atpar, tx_idx, nnofeat) {
  if (length(prior) > 1) {
    featuresnames <- names(atpar@features)[tx_idx[1:(length(tx_idx)-nnofeat)]]
    prior <- prior[featuresnames]
    prior <- c(prior, rep(0, nnofeat)) # setting prior of 0 to __no_feature
  } else if (length(prior) == 1) {
    prior <- rep(prior, length(tx_idx))
  }
  
  if (any(is.na(prior))) {
    stop(sprintf("Prior not found for feature: %s", 
                 featuresnames[which(is.na(prior))]))
  }
  return(prior)
}


## private function .atEstep()
## E-step of the EM algorithm based on Telescope
.atEstep <- function(Q, Theta, maskmulti, Pi) {
  # X <- t(t(Q) * Pi)
  # # X[maskmulti, ] <- t(t(X[maskmulti, ]) * Theta)
  # # quicker computation of previous line
  # wh <- which(X*maskmulti > 0, arr.ind = TRUE)
  # X[wh] <- t(t(X[wh]) * Theta[wh[, "col"]])
  X <- Q
  j <- rep(seq_len(ncol(X)), diff(X@p))
  X@x <- X@x * Pi[j]
  #wh <- which(X@x * maskmulti[X@i + 1] > 0)
  whmulti <- as.integer(maskmulti[X@i + 1])
  wh <- which(X@x * whmulti > 0)
  j <- rep(seq_len(ncol(X)), diff(X@p))
  X@x[wh] <- X@x[wh] * Theta[j][wh]
  X <- X[rowSums2(X) > 0, , drop = FALSE]
  X <- X/rowSums2(X)
  X
}

## private function .atMstepPi()
## M-step of the EM algorithm based on Telescope
.atMstepPi <- function(X, a) {
  Pi <- (colSums2(X) + a)/(sum(X@x) + sum(a))
  Pi
}

## private function .atMstepTheta()
## Update the estimate of the MAP value of θ
.atMstepTheta <- function(X, maskmulti, b) {
  Theta <- ((colSums2(X[maskmulti, , drop=FALSE]) + b) / 
              (sum(maskmulti) + sum(b)))
}

## private function .atFixedPointFun()
## fixed point function of the EM algorithm based on Telescope
.atFixedPointFun <- function(Pi, Thetaenv, Q, maskmulti, a, b) {
  Theta <- get("Theta", envir=Thetaenv)
  X <- .atEstep(Q, Theta, maskmulti, Pi)
  Pi2 <- .atMstepPi(X, a)
  Theta2 <- .atMstepTheta (X, maskmulti, b)
  assign("Theta", Theta2, envir=Thetaenv)
  
  Pi2
}
