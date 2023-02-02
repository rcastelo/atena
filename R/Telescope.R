#' Build a Telescope parameter object
#'
#' Build an object of the class \code{TelescopeParam}.
#'
#' @param bfl A \code{BamFile} or \code{BamFileList} object, or a character
#' string vector of BAM filenames.
#'
#' @param teFeatures A \code{GRanges} or \code{GRangesList} object. Elements
#' in this object should have names, which are used as a grouping factor for
#' genomic ranges forming a common locus (equivalent to "locus" column in
#' Telescope). This grouping is performed previous to TE expression
#' quantification, unlike the aggregation of quantifications performed when
#' the \code{aggregateby} parameter is specified, which is performed after
#' individual TE instances are quantified.
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
#' as in the original Telescope method: the overlap with the longest
#' overlapping length is kept.
#'
#' @param geneFeatures A \code{GRanges} or \code{GRangesList} object with the
#' gene annotated features to be quantified. The TEtranscripts approach for
#' gene expression quantification is used, in which overlaps with multi-mapping
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
#' @param fragments (Default FALSE) A logical; applied to paired-end data only.
#' When \code{fragments=FALSE}, the read-counting method only counts
#' ‘mated pairs’ from opposite strands (non-ambiguous properly paired reads), 
#' while when \code{fragments=TRUE} same-strand pairs, singletons, reads with 
#' unmapped pairs and other ambiguous or not properly paired fragments
#' are also counted (see "Pairing criteria" in 
#' \code{\link[GenomicAlignments]{readGAlignments}()}). \code{fragments=TRUE} 
#' is equivalent to the original Telescope algorithm. For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}()}.
#'
#' @param minOverlFract (Default 0.2) A numeric scalar. \code{minOverlFract}
#' is multiplied by the read length and the resulting value is used to
#' discard alignments for which the overlapping length (number of base
#' pairs the alignment and the feature overlap) is lower. When no minimum 
#' overlap is required, set \code{minOverlFract = 0}.
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
                            ovMode="ovUnion",
                            geneFeatures=NA,
                            singleEnd=TRUE,
                            strandMode=1L,
                            ignoreStrand=FALSE,
                            fragments=FALSE,
                            minOverlFract=0.2,
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
    
    new("TelescopeParam", bfl=bfl, features=features,
        aggregateby=aggregateby, ovMode=ovMode,
        singleEnd=singleEnd,ignoreStrand=ignoreStrand,
        strandMode=as.integer(strandMode), fragments=fragments,
        minOverlFract=minOverlFract, pi_prior=pi_prior,
        theta_prior=theta_prior, em_epsilon=em_epsilon,
        maxIter=as.integer(maxIter), reassign_mode=reassign_mode,
        conf_prob=conf_prob)
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
            
            cnt <- bplapply(x@bfl, .qtex_telescope, tspar=x, mode=x@ovMode,
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
#' @importFrom Matrix Matrix t which
#' @importFrom sparseMatrixStats rowSums2 colSums2
#' @importFrom SQUAREM squarem
.qtex_telescope <- function(bf, tspar, mode, yieldSize=1e6L) {
    mode <- match.fun(mode)
    readfun <- .getReadFunction(tspar@singleEnd, tspar@fragments)
    .checkreassignModes(tspar)
    sbflags <- .getScanBamFlag_ts(tspar@fragments)
    param <- ScanBamParam(flag=sbflags, what="flag", tag="AS")
    iste <- as.vector(attributes(tspar@features)$isTE[,1])
    if (any(duplicated(names(tspar@features[iste])))) {
        stop(".qtex_telescope: transposable element annotations do not contain unique names for each element")
    }
    ov <- Hits(nLnode=0, nRnode=length(tspar@features), sort.by.query=TRUE)
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
                                list(strandMode=tspar@strandMode)[strand_arg],
                                list(use.names=TRUE))))) {
        alnreadids <- c(alnreadids, names(alnreads))
        asvalues <- c(asvalues, .getAlignmentASScoreTS(alnreads, tag = "AS"))
        readlen <- .getAlignmentLength(alnreads)
        alen <- c(alen, readlen)
        salnmask <- c(salnmask, any(.secondaryAlignmentMask(alnreads)))
        thisov <- mode(alnreads, tspar@features,
                       ignoreStrand=tspar@ignoreStrand, inter.feature=FALSE)
        
        # Selecting the overlaps with min overlap higher than threshold
        ovlength <- .getOverlapLength(alnreads, thisov, tspar)
        yesminov <- ovlength > tspar@minOverlFract*readlen[queryHits(thisov)]
        thisov <- thisov[yesminov]
        ovlength <- ovlength[yesminov]
        
        # Selecting the best overlap for alignments overlapping > 1 feature
        multiov <- (duplicated(queryHits(thisov)) |
                        duplicated(queryHits(thisov), fromLast = TRUE))
        if (any(multiov)) {
          int <- ovlength[multiov]
          intmax <- aggregate(int, by = list(queryHits(thisov)[multiov]),
                              FUN = which.max)
          whpos <- which(!duplicated(queryHits(thisov)[multiov]))
          whpos <- whpos - 1
          intmax <- intmax[!duplicated(intmax$Group.1),]
          thisov <- thisov[-which(multiov)[-(whpos + intmax$x)]]
        }

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
    cntvec <- .tsEMstep(tspar, alnreadids, readids, ov, asvalues,
                        iste, maskuniqaln, mt, alen)
    setNames(as.integer(cntvec), names(cntvec))
}

## private function .getNoFeatureOv(), obtains the overlaps to "no_feature"
#' @importFrom S4Vectors Hits queryHits subjectHits nRnode nLnode from to
.getNoFeatureOv <- function(maskuniqaln, ov, alnreadids) {
    if (!all(maskuniqaln)) {
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
    } else {
      ovnofeat <- ov
    }
    return(ovnofeat)
}

#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom Matrix Matrix t which
#' @importFrom SQUAREM squarem
#' @importClassesFrom Matrix lgCMatrix
.tsEMstep <- function(tspar, alnreadids, readids, ov, asvalues,
                      iste, maskuniqaln, mt, alen) {
    ## initialize vector of counts derived from multi-mapping reads
    cntvec <- rep(0L, length(tspar@features) + 1)
    
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
    asvalues <- .rescaleAS(asvalues, alen = alen)
    QmatTS <- .buildOvValuesMatrix(tspar, ov, asvalues, alnreadidx,
                                    rd_idx, tx_idx)
    
    if (!all(iste)) { ## Correcting for preference of reads to genes/TEs
        if (all(maskuniqaln)) {
          QmatTS <- .correctPreferenceTS(QmatTS, maskuniqaln, mt, istex)
          stopifnot(!any(rowSums2(QmatTS[,istex]) > 0 &
                           rowSums2(QmatTS[,!istex]) > 0))
        } else {
            QmatTS_nofeat <- QmatTS[,ncol(QmatTS), drop = FALSE]
            QmatTS <- .correctPreferenceTS(QmatTS[,-ncol(QmatTS)], maskuniqaln,
                                           mt, istex)
            stopifnot(!any(rowSums2(QmatTS[,istex]) > 0 &
                               rowSums2(QmatTS[,!istex]) > 0))
            QmatTS <- cbind(QmatTS, QmatTS_nofeat)
        }
    }
    
    
    # --- EM-step --- 
    # QmatTS <- QmatTS / rowSums2(QmatTS)
    QmatTS@x <- QmatTS@x / rowSums2(QmatTS)[QmatTS@i +1]
    
    # Getting 'y' indicator of unique and multi-mapping status from QmatTS 
    # maskmulti <- ifelse(rowSums2(QmatTS > 0) == 1, 0, 1)
    maskmulti <- rowSums2(QmatTS > 0) > 1
    
    # Telescope (Bendall et al.(2019)) defines the initial π estimate uniformly
    PiTS <- rep(1 / length(tx_idx), length(tx_idx))
    
    # Telescope defines an additional reassignment parameter (θ) as uniform
    Theta <- rep(1 / length(tx_idx), length(tx_idx))
    
    ## The SQUAREM algorithm to run the EM procedure
    a <- as.numeric(tspar@pi_prior) # 0
    b <- as.numeric(tspar@theta_prior) # 0
    Thetaenv <- new.env()
    assign("Theta", Theta, envir=Thetaenv)
    tsres <- squarem(par=PiTS, Thetaenv=Thetaenv, Q=QmatTS,maskmulti=maskmulti,
                    a=a, b=b, fixptfn=.tsFixedPointFun,
                    control=list(tol=tspar@em_epsilon, maxiter=tspar@maxIter))
    PiTS <- tsres$par
    PiTS[PiTS < 0] <- 0 ## Pi estimates are sometimes negatively close to zero
    # --- end EM-step ---
    
    Theta <- get("Theta", envir=Thetaenv)
    X <- .tsEstep(QmatTS, Theta, maskmulti, PiTS)
    cntvec <- .reassign(X, tspar@reassign_mode, tspar@conf_prob, cntvec, tx_idx)
    
    names(cntvec) <- c(names(tspar@features), "no_feature")
    nofeat <- cntvec["no_feature"]
    cntvec <- .tssummarizeCounts(cntvec[-length(cntvec)], iste, tspar)
    cntvec <- c(cntvec, nofeat)
    cntvec
}



## private function .checkreassignModes()
## Checks that 'reassign_mode' is valid
.checkreassignModes <- function(tspar) {
  if (!(tspar@reassign_mode %in% c("exclude","choose","average","conf")))
    stop("'reassign_mode' must be one of: 'exclude', 'choose', 'average' or 'conf'")
}


## private function .reassign()
## Implements 4 different reassigning methods from Telescope (exclude, choose,
## average and conf)
#' @importFrom sparseMatrixStats rowSums2 colSums2
#' @importFrom Matrix summary
#' @importFrom stats runif
.reassign <- function(X, reassign_mode, conf_prob, cntvec, tx_idx) {
  
  if (reassign_mode %in% c("exclude", "choose", "average")) {
    maxbyrow <- rowMaxs(X)
    # Older versions
    # Xind <- X == maxbyrow
    # Xind <- (X / maxbyrow) == 1
    Xind <- as(as(as(X, "lMatrix"), "generalMatrix"), "CsparseMatrix")
    Xind@x <- (X@x /maxbyrow[X@i+1]) == 1
    nmaxbyrow <- rowSums2(Xind)
    cntvec[tx_idx] <- colSums2(Xind[nmaxbyrow == 1, ])
  
    if (reassign_mode == "choose" & any(nmaxbyrow > 1)) {
      Xind2 <- Xind[nmaxbyrow > 1, ]
      Xind2_s <- summary(Xind2)
      Xind2_s <- Xind2_s[Xind2_s$x == TRUE,]
      Xind2_s <- Xind2_s[order(Xind2_s[,"i"]),]
      bestov <- as.vector(table(Xind2_s[,"i"]))
      # selected_ov <- vapply(bestov, FUN = function(x) sample(x = x, size = 1), 
      #                       FUN.VALUE = integer(1L))
      selected_ov <- ceiling(runif(n=nrow(Xind2))*rowSums2(Xind2))
      ovcol <- table(Xind2_s[c(0,cumsum(bestov)[-length(bestov)]) + selected_ov,"j"])
      whcol <- as.integer(names(ovcol))
      cntvec[tx_idx][whcol] <- cntvec[tx_idx][whcol] + as.integer(ovcol)
    }
    
    if (reassign_mode == "average" & any(nmaxbyrow > 1)) {
      Xind2 <- Xind[nmaxbyrow > 1, ]
      cntvec[tx_idx] <- cntvec[tx_idx] + colSums2(Xind2/rowSums2(Xind2))
    }
    
  } else if (reassign_mode == "conf") {
    whbelow <- X@x < conf_prob
    X@x[whbelow] <- 0
    X_conf <- X[rowSums2(X) > 0,]
    # X_conf@i indicates row and is 0-based
    # Faster sparse implementation of X_conf/rowSums2(X_conf)
    X_conf@x <- X_conf@x / rowSums2(X_conf)[X_conf@i +1] 
    cntvec[tx_idx] <- colSums2(X_conf)
    
  } else {
    stop("'reassign_mode' should be one of 'exclude', 'choose', 'average' or 'conf'")
  }
  
  cntvec
}


#' @importFrom scales rescale
.rescaleAS <- function(asvalues, alen) {
  asrescale <- rescale(asvalues, c(1, max(asvalues) - min(asvalues) + 1))
  asrescale <- asrescale + alen
  asrescale <- asrescale/max(asrescale)
  asrescale <- expm1(asrescale*100)
  asrescale
}


## private function .correctPreferenceTS()
## Corrects QmatTS for preference of unique/multi-mapping reads to
## genes/TEs, respectively, in Telescope
.correctPreferenceTS <- function(QmatTS, maskuniqaln, mt, istex) {
    indx <- (rowSums2(QmatTS[,istex]) > 0) & (rowSums2(QmatTS[,!istex]) > 0)
    
    ## Assigning unique reads mapping to both a TE and a gene as gene counts
    # Which unique reads overlap to both genes and TEs?
    idxu <- indx & maskuniqaln[mt]
    if (any(idxu)) {
        # QmatTS[idxu,istex] <- FALSE
        whu <- which(QmatTS[idxu,istex] > 0, arr.ind = TRUE)
        if (is(whu, "integer")) {
            whu <- as.matrix(data.frame(row = 1, col = whu))
        }
        whudf <- cbind(which(idxu)[whu[,"row"]], which(istex)[whu[,"col"]])
        QmatTS[whudf] <- FALSE
    }
    
    ## Removing overlaps of multi-mapping reads to genes if at least one
    ## alignment of the read overlaps a TE
    idxm <- indx & !maskuniqaln[mt]
    if (any(idxm)) {
        # QmatTS[idxm,!istex] <- FALSE
        whm <- which(QmatTS[idxm,!istex] > 0, arr.ind = TRUE)
        if (is(whm, "integer")) {
          whm <- as.matrix(data.frame(row = 1, col = whm))
        }
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


## private function .tsMstepPi()
## M-step of the EM algorithm of Telescope
.tsMstepPi <- function(X, a) {
  Pi <- (colSums2(X) + a)/(sum(X@x) + a * ncol(X))
    Pi
}

## private function .tsMstepTheta()
## Update the estimate of the MAP value of θ
.tsMstepTheta <- function(X, maskmulti, b) {
    Theta <- ((colSums2(X[maskmulti, , drop=FALSE]) + b) / 
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


#' @importFrom S4Vectors mcols first second
.getAlignmentASScoreTS <- function(aln, tag) {
  if (is(aln, "GAlignments"))
    score <- as.integer(mcols(aln)[[tag]])
  else if (is(aln, "GAlignmentPairs")) {
    ## take the sum of the score of each mate
    score <- as.integer(mcols(first(aln))[[tag]]) + 
        as.integer(mcols(second(aln))[[tag]])
  } else if (is(aln, "GAlignmentsList")) {
    l <- lengths(aln)
    score <- aggregate(mcols(unlist(aln, use.names = FALSE))[[tag]],
                        by = list(rep(seq_along(l),l)), FUN = sum)
    score <- score$x
    mate_status <- mcols(aln)$mate_status == "mated"
    score[!mate_status] <- unlist(lapply(aln[!mate_status], 
                                  function(x) mean(mcols(x)[,tag])*2))
                                  # multiplied by 2 since there are 2 mates
  } else {
    stop(sprintf(".getAlignmentTagScore: wrong class %s\n", class(aln)))
  }
  as.integer(score)
}

#' @importFrom S4Vectors first second
#' @importFrom GenomicRanges granges
#' @importFrom GenomicAlignments seqnames
.getAlignmentLength <- function(alnreads) {
  if (is(alnreads, "GAlignments"))
    readlen <- qwidth(alnreads)
  else if (is(alnreads, "GAlignmentPairs")) {
      ## take the length of the region comprised by the 2 mates
      readlen <- width(granges(alnreads))
  } else if (is(alnreads, "GAlignmentsList")) {
      # readlen <- width(granges(alnreads, ignore.strand=TRUE))
      # # In case of reads with space between the two mates, the read length
      # # assigned corresponds to twice the maximum alignment length, to account
      # # for the length of the two mates, but not the region between them
      # maxlen <- max(qwidth(unlist(alnreads)))
      # readlen[readlen>maxlen] <- maxlen*2
      
      lsum <- sum(width(alnreads))
      # unmated alignments
      unmat <- mcols(alnreads)$mate_status == "unmated"
      # pairs with discordant chr
      unmat <- unmat | lengths(unique(seqnames(alnreads)))>1
      # Unmated alignments the assigned read length is the median read length
      lsum[unmat] <- median(width(unlist(alnreads)))
      ltogether <- integer(length = length(alnreads))
      ltogether[!unmat] <- width(granges(alnreads[!unmat], ignore.strand=TRUE))
      ltogether[unmat] <- median(width(unlist(alnreads)))
      ltogether[ltogether > lsum] <- lsum[ltogether > lsum]
      readlen <- ltogether
      
  } else {
    stop(sprintf(".getAlignmentLength: wrong class %s\n", class(alnreads)))
  }
  readlen
}


#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom GenomicRanges pintersect
#' @importFrom GenomeInfoDb seqlevels<- seqlevels
.getOverlapLength <- function(alnreads, thisov, tspar) {
    features <- tspar@features
    seqlev <- unique(c(seqlevels(features), seqlevels(alnreads)))
    seqlevels(features) <- seqlev
    seqlevels(alnreads) <- seqlev
    
    if (is(alnreads, "GAlignmentsList")) {
        l <- lengths(alnreads)[queryHits(thisov)]
        features_ov <- rep(features[subjectHits(thisov)], l)
        ovlength <- width(pintersect(GRanges(unlist(alnreads[queryHits(thisov)])),
                                      features_ov,
                                      ignore.strand = tspar@ignoreStrand,
                                      strict.strand=FALSE))
        if (is(ovlength, "CompressedIntegerList")) {
          ovlength <- max(ovlength)
        }
        ovlength_ag <- aggregate(ovlength, by = list(rep(seq_along(l), l)), 
                                  FUN = max)
        ovlength <- ovlength_ag$x
    } else {
        ovlength <- width(pintersect(GRanges(alnreads[queryHits(thisov)]),
                                      features[subjectHits(thisov)],
                                      ignore.strand = tspar@ignoreStrand,
                                      strict.strand=FALSE))
        if (is(ovlength, "CompressedIntegerList")) {
          ovlength <- max(ovlength)
        }
    }
    ovlength
}

#' @importFrom Rsamtools ScanBamParam
.getScanBamFlag_ts <- function(fragments) {
  if (fragments == TRUE) {
    sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                           isDuplicate=FALSE,
                           isNotPassingQualityControls=FALSE)
  } else {
    sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                           isDuplicate=FALSE,
                           isNotPassingQualityControls=FALSE,
                           isProperPair=TRUE)
  }
  sbflags
}
