library(Rsamtools)
library(Matrix)
library(sparseMatrixStats)

test_new_Telescope_functions <- function(){
  ####################################################################
  ######### Creating all params needed to check functions ############
  ####################################################################
  
  bamfiles <- list.files(system.file("extdata", package="atena"),
                         pattern="*.bam", full.names=TRUE)
  TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                         package="atena"))
  tspar <- TelescopeParam(bamfiles, teFeatures = TE_annot,
                          singleEnd = TRUE, ignoreStrand=TRUE,
                          aggregateby = c("repName"))        
  mode <- ovUnion
  readfun <- atena:::.getReadFunction(tspar@singleEnd, tspar@fragments)
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
  bf <- tspar@bfl[[1]]
  yieldSize <- 1e6L
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
    ov <- atena:::.appendHits(ov, thisov)
  }
  on.exit(close(bf))
  maskuniqaln <- !(duplicated(alnreadids) | 
                     duplicated(alnreadids, fromLast = TRUE))
  readids <- unique(alnreadids[queryHits(ov)])
  ov <- atena:::.getNoFeatureOv(maskuniqaln, ov, alnreadids)
  mt <- match(readids, alnreadids)
  readids <- unique(alnreadids[queryHits(ov)]) # updating 'readids'
  cntvec <- rep(0L, length(tspar@features) + 1)
  alnreadidx <- match(alnreadids, readids)
  rd_idx <- sort(unique(alnreadidx[queryHits(ov)]))
  
  tx_idx <- sort(unique(subjectHits(ov)))
  istex <- as.vector(iste[tx_idx])[-length(tx_idx)]  # removing "no_feature"
  
  asvalues <- (asvalues-min(asvalues)+1) / (max(asvalues)+1 - min(asvalues))
  QmatTS <- atena:::.buildOvValuesMatrix(tspar, ov, asvalues, alnreadidx,
                                         rd_idx, tx_idx)
  QmatTS@x <- QmatTS@x / rowSums2(QmatTS)[QmatTS@i +1]
  maskmulti <- ifelse(rowSums2(QmatTS > 0) == 1, 0, 1)
  PiTS <- rep(1 / length(tx_idx), length(tx_idx))
  Theta <- rep(1 / length(tx_idx), length(tx_idx))
  
  ####################################################################
  ############################  TEST 1 ###############################
  ####################################################################
  
  
  old_tsEstep <- function() {
    X <- t(t(QmatTS) * PiTS)
    wh <- which(X*maskmulti > 0, arr.ind = TRUE)
    X[wh] <- t(t(X[wh]) * Theta[wh[, "col"]])
    X <- X[rowSums(X)>0,, drop=FALSE]
    X <- X / rowSums(X)
    X
  }
  
  new_tsEstep <- function() {
    X <- QmatTS
    j <- rep(1:ncol(X), diff(X@p))
    X@x <- X@x * PiTS[j]
    wh <- which(X@x * maskmulti[X@i + 1] > 0)
    X@x[wh] <- X@x[wh] * Theta[j][wh]
    X <- X[rowSums2(X) > 0, , drop = FALSE]
    X <- X/rowSums2(X)
    X
  }
  
  checkEquals(
    old_tsEstep(),
    new_tsEstep()
  )
  
  ####################################################################
  ############################  TEST 2 ###############################
  ####################################################################
  
  X <- new_tsEstep()
  maxbyrow <- rowMaxs(X)
  
  old_Xind <- function(){
    Xind <- (X / maxbyrow) == 1
  }
  
  new_Xind <- function(){
    Xind <- as(as(as(X, "lMatrix"), "generalMatrix"), "CsparseMatrix")
    Xind@x <- (X@x /maxbyrow[X@i+1]) == 1
    Xind
  }
  
  checkEquals(
    old_Xind(),
    new_Xind()
  )
}
  
