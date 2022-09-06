library(Rsamtools)
library(Matrix)
library(sparseMatrixStats)

test_new_TEtranscript_functions <- function(){
  
  ####################################################################
  ######### Creating all params needed to check functions ############
  ####################################################################
  
  bamfiles <- list.files(system.file("extdata", package="atena"),
                         pattern="*.bam", full.names=TRUE)
  TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                         package="atena"))
  ttpar <- TEtranscriptsParam(bamfiles, teFeatures = TE_annot,
                              singleEnd = TRUE, ignoreStrand=TRUE,
                              aggregateby = c("repName"))        
  mode=match.fun(ovUnion)
  readfun <- atena:::.getReadFunction(ttpar@singleEnd, ttpar@fragments)
  sbflags <- scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE,
                         isNotPassingQualityControls=FALSE)
  param <- ScanBamParam(flag=sbflags, tag="AS")
  iste <- as.vector(attributes(ttpar@features)$isTE[,1])
  if (any(duplicated(names(ttpar@features[iste])))) {
    stop(".qtex_tetranscripts: transposable element annotations
         do not contain unique names for each element")
  }
  ov <- Hits(nLnode=0, nRnode=length(ttpar@features), sort.by.query=TRUE)
  alnreadids <- character(0)
  avgreadlen <- integer(0)
  d <- numeric(0)
  strand_arg <- "strandMode" %in% formalArgs(readfun)
  yieldSize <- 1e6L
  bf <- ttpar@bfl[[1]]
  yieldSize(bf) <- yieldSize
  open(bf)
  while (length(alnreads <- do.call(readfun, 
                                    c(list(file = bf), list(param=param),
                                      list(strandMode=ttpar@strandMode)[strand_arg],
                                      list(use.names=TRUE))))) {
    avgreadlen <- c(avgreadlen, atena:::.getAveLen(ttpar, alnreads))
    if (ttpar@singleEnd == FALSE) {
      d <- c(d, abs(start(first(alnreads)) - start(second(alnreads))))
    }
    alnreadids <- c(alnreadids, names(alnreads))
    thisov <- mode(alnreads, ttpar@features, minOverlFract=0L,
                   ignoreStrand=ttpar@ignoreStrand)
    ov <- atena:::.appendHits(ov, thisov)
  }
  on.exit(close(bf))
  
  maskuniqaln <- !(duplicated(alnreadids) | duplicated(alnreadids, 
                                                       fromLast = TRUE))
  readids <- unique(alnreadids[queryHits(ov)])
  tx_idx <- sort(unique(subjectHits(ov)))
  
  ####################################################################
  ############################  TEST 1 ###############################
  ####################################################################
  
  old_buildOvAlignmentsMatrix <- function(ov, arids, rids, fidx){
    oamat <- Matrix(FALSE, nrow = length(rids), ncol = length(fidx))
    mt1 <- match(arids[queryHits(ov)], rids)
    mt2 <- match(subjectHits(ov), fidx)
    oamat[cbind(mt1, mt2)] <- TRUE
    oamat
  }
  
  new_buildOvAlignmentsMatrix <-  function(ov, arids, rids, fidx){
    mt1 <- match(arids[queryHits(ov)], rids)
    mt2 <- match(subjectHits(ov), fidx)
    positions <- cbind(mt1, mt2)
    oamat <- sparseMatrix(positions[,1], positions[,2], x=TRUE)
    oamat
  }
  
  checkEquals(old_buildOvAlignmentsMatrix(ov, alnreadids, readids, tx_idx),
              new_buildOvAlignmentsMatrix(ov, alnreadids, readids, tx_idx))
  
  ####################################################################
  ############################  TEST 2 ###############################
  ####################################################################
  
  ovalnmat <- atena:::.buildOvAlignmentsMatrix(ov, alnreadids, readids, tx_idx)
  mt <- match(readids, alnreadids)
  istex <- as.vector(iste[tx_idx])
  ovalnmat <- ovalnmat[!maskuniqaln[mt], istex]
  yesov <- rowSums2(ovalnmat)>0
  ovalnmat <- ovalnmat[yesov,]
  readids <- readids[!maskuniqaln[mt]][yesov]
  
  old_Qmat_generation <- function (){
    Qmat <- Matrix(0, nrow=length(readids), ncol=length(tx_idx[istex]),
                   dimnames=list(readids, NULL))
    Qmat[which(ovalnmat, arr.ind=TRUE)] <- 1
    Qmat
  }
  
  new_Qmat_generation <- function(){
    x <- as.integer(ovalnmat@x)
    Qmat <- sparseMatrix(i=ovalnmat@i, p=ovalnmat@p, x=x, index1 = FALSE,
                         dimnames =  list(readids, NULL))
  }
  
  checkEquals(old_Qmat_generation(), new_Qmat_generation(), 
              msg = "Testing Qmat generation")
  
  ####################################################################
  ############################  TEST 3 ###############################
  ####################################################################
  
  Qmat <- new_Qmat_generation()
  
  old_divide_Qmat <- function(){
    Qmat / rowSums2(ovalnmat)
  }
  
  new_divide_Qmat <- function(){
    qmat <- Qmat
    qmat@x <- Qmat@x / rowSums2(ovalnmat)[Qmat@i +1]
    qmat
  }
  
  checkEquals(old_divide_Qmat(), new_divide_Qmat(),
              msg = "Testing dividing Qmat by its rowSums")
  
  ####################################################################
  ############################  TEST 4 ###############################
  ####################################################################
  
  Pi <- colSums2(Qmat)
  Pi <- Pi / sum(Pi)
  probmassbyread <- as.vector(ovalnmat %*% Pi)
  
  old_cntvecovtx <- function(){
    wh <- which(ovalnmat, arr.ind=TRUE)
    cntvecovtx <- rep(0, ncol(ovalnmat)) #cntvecovtx <- rep(0, length(tx_idx))
    x <- tapply(Pi[wh[, "col"]] / probmassbyread[wh[, "row"]], wh[, "col"],
                FUN=sum, na.rm=TRUE)
    cntvecovtx[as.integer(names(x))] <- x
    cntvecovtx
  }
  
  new_cntvecovtx <- function(){
    ovalnmat2 <- as(as(as(ovalnmat, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    i <- ovalnmat2@i + 1
    ovalnmat2@x <- ovalnmat2@x / probmassbyread[i]
    j <- rep(1:ncol(ovalnmat2), diff(ovalnmat2@p))
    ovalnmat2@x <- ovalnmat2@x * Pi[j]
    cntvecovtx <- sparseMatrixStats::colSums2(ovalnmat2, na.rm=TRUE)
    cntvecovtx
  }
  
  checkEquals(old_cntvecovtx(), new_cntvecovtx(), 
              msg = "Testing generating 'cntvecovtx' vector")
  
  ####################################################################
  ############################  TEST 4 ###############################
  ####################################################################
  
  Q <- rsparsematrix(1000, 200, 0.2)
  Pi <- colSums2(Q)
  
  old_tsEstep <- function(){
    X <- t(t(Q) * Pi)
    X <- X[rowSums(X) > 0, , drop = FALSE]
    X <- X/rowSums(X)
    X
  }
  
  new_tsEstep <- function(){
    X <- Q
    j <- rep(1:ncol(X), diff(X@p))
    X@x <- X@x * Pi[j]
    X <- X[rowSums2(X)>0,, drop=FALSE]
    X <- X / rowSums2(X)
    X
  }
  
  checkEquals(old_tsEstep(), new_tsEstep(), 
              msg = "Testing .tsEstep function")
  
}
