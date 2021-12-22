.h5_ttQuantExpress <- function(ov, alnreadids, readids, tx_idx, ttpar, iste,
                            maskuniqaln, avgreadlen) {
  
  ## build a matrix representation of the overlapping alignments
  ovalnh5 <- .buildOvAlignmentsH5(ov, alnreadids, readids, tx_idx)
  
  mt <- match(readids, alnreadids)
  multigcnt <- rep(0L, length(ttpar@features))
  istex <- as.vector(iste[tx_idx])
  
  if (!all(iste)) {
    ## Correcting for preference of unique/multi-mapping reads to genes or TEs
    ovalnh5 <- .h5_correctPreference(ovalnh5, maskuniqaln, mt, istex)
  }
  
  # Getting counts from unique reads
  uniqcnt <- .h5_countUniqueRead(ttpar, ovalnh5, maskuniqaln, mt, tx_idx, istex)
  
  if (!all(iste) & any(!maskuniqaln)) {
    ## Getting gene counts where a multimapping read maps to > 1 gene
    multigcnt <- .countMultiReadsGenes(ttpar, ovalnh5, maskuniqaln, mt,
                                       iste, istex, tx_idx, readids, 
                                       alnreadids, ov, uniqcnt)
  }
  cntvec <- rep(0L, length(ttpar@features))
  cntvec <- .h5_ttEMstep(maskuniqaln, mt, ovalnh5, istex, tx_idx, readids, ttpar,
                      avgreadlen, cntvec)
  ## Summarize counts of unique and multimapping reads
  cntvec <- .summarizeCounts(iste, cntvec, uniqcnt, multigcnt, ttpar)
  cntvec
}


#' @importFrom Matrix which
#' @importFrom SQUAREM squarem
#' @importFrom IRanges ranges
#' @importFrom DelayedArray t
#' @importFrom DelayedMatrixStats rowSums2 colSums2
#' @importFrom HDF5Array h5writeDimnames
.h5_ttEMstep <- function(maskuniqaln, mt, ovalnh5, istex, tx_idx, readids, ttpar,
                      avgreadlen, cntvec) {
  if (sum(!maskuniqaln[mt]) > 0) { ## multi-mapping reads
    ## TEtranscripts doesn't use uniquely-aligned reads to inform the
    ## procedure of distributing multiple-mapping reads, as explained in
    ## Jin et al. (2015) pg. 3594, "to reduce potential bias to certain TEs."
    ## Hence, once counted, we discard unique alignments
    ovalnh5 <- ovalnh5[!maskuniqaln[mt], istex]
    yesov <- rowSums2(ovalnh5)>0
    ovalnh5 <- ovalnh5[yesov,]
    readids <- readids[!maskuniqaln[mt]][yesov]
    ## the Qmat_h5 matrix stores row-wise the probability that read i maps to
    ## a transcript j, assume uniform probabilities by now
    print("Creating an hdf5 Qmat_h5 matrix")
    Qmat_h5_file <- tempfile(fileext = ".h5")
    h5createFile(Qmat_h5_file)
    h5createDataset(Qmat_h5_file, dataset = "data", dims = c(length(readids), length(tx_idx)),
                    chunk = c(100L, 100L), fillValue = 0, level = 1)
    h5writeDimnames(list(readids, NULL), Qmat_h5_file, "data")
    Qmat_h5 <- HDF5Array(Qmat_h5_file, name="data")
    
    ## update Qmat_h5_h5
    ma2 <- which(ovalnh5, arr.ind = TRUE)
    updateElements(Qmat_h5_file, "data", ma2, 1, TRUE)
    
    Qmat_h5 <- Qmat_h5 / rowSums2(ovalnh5)
    
    ## Pi, corresponding to rho in Equations (1), (2) and (3) in Jin et al.
    ## (2015) stores probabilities of expression for each transcript, corrected
    ## for its effective length as defined in Eq. (1) of Jin et al. (2015)
    Pi <- colSums2(Qmat_h5)
    
    if (is(ttpar@features,"GRangesList")) {
      elen <- as.numeric(width(ttpar@features[tx_idx][istex])) - avgreadlen+1
    } else {
      elen <- width(ttpar@features[tx_idx][istex]) - avgreadlen + 1
    }
    Pi <- .correctForTxEffectiveLength(Pi, elen)
    
    ## use the SQUAREM algorithm to achieve faster EM convergence
    emres <- squarem(par=Pi, Q=Qmat_h5, elen=elen,
                     fixptfn=.h5_ttFixedPointFun,
                     control=list(tol=ttpar@tolerance, maxiter=ttpar@maxIter))
    Pi <- emres$par
    Pi[Pi < 0] <- 0 ## Pi estimates are sometimes negatively close to zero
    Pi <- Pi / sum(Pi)
    ## use the estimated transcript expression probabilities
    ## to finally distribute ambiguously mapping reads
    probmassbyread <- as.vector(ovalnh5 %*% Pi) 
    # cntvecovtx <- rowSums(t(ovalnh5 / probmassbyread) * Pi, na.rm=TRUE)
    # sparce version of the previous commented line, improves memory usage
    wh <- which(ovalnh5, arr.ind=TRUE)
    cntvecovtx <- rep(0, ncol(ovalnh5)) #cntvecovtx <- rep(0, length(tx_idx))
    x <- tapply(Pi[wh[, "col"]] / probmassbyread[wh[, "row"]], wh[, "col"],
                FUN=sum, na.rm=TRUE)
    cntvecovtx[as.integer(names(x))] <- x
    cntvec[tx_idx][istex] <- cntvecovtx
  }
  cntvec
}


## function to update an h5 using a matrix of positions 
## solution provided by Mike Smith in this slack conversation:
## https://tinyurl.com/37c4wnsx
#' @importFrom rhdf5 H5Fopen H5Fclose H5Dopen H5Dclose H5Dget_space H5Sclose H5Screate_simple H5Dwrite
updateElements <- function(h5_file, dataset_name, row_col_matrix, value, isValueRep = TRUE) {
  
  ## Open HDF5 file, dataset,  dataspace
  ## Important they're closed when we exit the function
  fid <- H5Fopen(h5_file)
  on.exit(H5Fclose(fid))
  did <- H5Dopen(fid, name = dataset_name)
  on.exit(H5Dclose(did), add = TRUE)
  file_dspace <- H5Dget_space(did)
  on.exit(H5Sclose(file_dspace), add = TRUE)
  
  ## There's a lot going on here.  We need to transform from row,col to col,row
  ## and make into a 1D stream of coordinates. See
  ## https://portal.hdfgroup.org/display/HDF5/H5S_SELECT_ELEMENTS for more info
  coords <- as.integer( rev(t(row_col_matrix)) )
  numElements <- nrow(row_col_matrix)
  
  ## Set selection on dataspace to be only the points we're interested in
  res <- .Call("_H5Sselect_elements", file_dspace@ID, 0L, numElements, coords, PACKAGE='rhdf5')
  if(res < 0) { stop("Selection failed") }
  
  ## create an HDF5 dataspace representing the data we're writing
  mem_dspace <- H5Screate_simple(dims = numElements, maxdims = NULL)
  on.exit(H5Sclose(mem_dspace), add = TRUE)
  
  ## perform writing
  if(isValueRep){
    H5Dwrite(h5dataset = did, buf = rep(value, numElements),
             h5spaceMem = mem_dspace, h5spaceFile = file_dspace) 
  } else {
    H5Dwrite(h5dataset = did, buf = value,
             h5spaceMem = mem_dspace, h5spaceFile = file_dspace)
  }
}

## build a matrix representation of the overlapping alignments
## filled with FALSE values in an h5 file and then updated with
## values coming from a position matrix
#' @importFrom rhdf5 h5createFile h5createDataset
#' @importFrom HDF5Array HDF5Array
.buildOvAlignmentsH5 <- function(ov, arids, rids, fidx){
  # create HDF5 file full of 'FALSE'
  my_hdf5_file <- tempfile(fileext = ".h5")
  h5createFile(my_hdf5_file)
  
  h5createDataset(my_hdf5_file, dataset = "data", dims = c(length(rids), length(fidx)),
                  chunk = c(100L, 100L), storage.mode = "logical", fillValue = FALSE)
  
  # create position matrix
  mt1 <- match(arids[queryHits(ov)], rids)
  mt2 <- match(subjectHits(ov), fidx)
  ma <- cbind(mt1, mt2)
  
  # update h5 matrix
  updateElements(my_hdf5_file, "data", ma, TRUE, TRUE)
  return(HDF5Array(my_hdf5_file, name="data"))
}

## private function .h5_ttEstep()
## E-step of the EM algorithm of TEtranscripts
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom DelayedArray read_block write_block rowAutoGrid
#' @importClassesFrom DelayedArray DelayedArray RealizationSink
#' @importClassesFrom Matrix sparseMatrix
.h5_ttEstep <- function (Q, Pi) {
  
  print("Calculating ttEstep()")
  sink <- HDF5RealizationSink(dim = dim(Q), 
                              chunkdim = c(100L, 100L), level=2)
  grid <- rowAutoGrid(Q, nrow = 1000L)
  
  for(bid in seq_along(grid)){
    viewport <- grid[[bid]]
    block <- read_block(Q, viewport, as.sparse=TRUE)
    block <- as(block, "sparseMatrix")
    j <- rep(1:ncol(block), diff(block@p))
    block@x <- block@x * Pi[j] 
    sink <- write_block(sink, viewport, block)
    message(bid, " block of ", max(seq_along(grid)), " blocks.")
  }
  close(sink)
  X <- as(sink, "DelayedArray")
  
  sumrow1 <- DelayedMatrixStats::rowSums2(X)
  X <- X[sumrow1 > 0, , drop = FALSE]
  
  sumrow2 <- DelayedMatrixStats::rowSums2(X)  
  X <- X/sumrow2
  X
  
  ## leave X as a DelayedArray, since next step (colSums2()) goes pretty fast
  ## (no need to re-write X as an HDF5Array)
}

## private function .h5_ttMstep()
## M-step of the EM algorithm of TEtranscripts
.h5_ttMstep <- function (X) {
  print("Calculating ttMstep()")
  sumCols <- colSums2(X)
  sumOfSumCols <- sum(sumCols)
  Pi <- sumCols/sumOfSumCols
  Pi
}

## private function .h5_ttFixedPointFun()
## fixed point function of the EM algorithm of TEtranscripts
.h5_ttFixedPointFun <- function (Pi, Q, elen) {
  X <- .h5_ttEstep(Q, Pi)
  Pi2 <- .h5_ttMstep(X) 
  Pi2 <- .correctForTxEffectiveLength(Pi2, elen)
  Pi2
}


## private function .h5_correctPreference()
## Corrects ovalnh5 for preference of unique/multi-mapping reads to
## genes/TEs, respectively
#' @importFrom DelayedArray cbind which
.h5_correctPreference <- function(ovalnh5, maskuniqaln, mt, istex) {
  indx <- (rowSums2(ovalnh5[,istex]) > 0) & (rowSums2(ovalnh5[,!istex]) > 0)
  
  ## Assigning unique reads mapping to both a TE and a gene as gene counts
  # Which unique reads overlap to both genes and TEs?
  idxu <- indx & maskuniqaln[mt]
  if (any(idxu)) {
    # ovalnh5[idxu,istex] <- FALSE
    whu <- which(ovalnh5[idxu,istex], arr.ind = TRUE)
    whudf <- cbind(which(idxu)[whu[,"row"]], which(istex)[whu[,"col"]])
    ovalnh5[whudf] <- FALSE
  }
  
  ## Removing overlaps of multi-mapping reads to genes if at least one
  ## alignment of the read overlaps a TE
  idxm <- indx & !maskuniqaln[mt]
  if (any(idxm)) {
    # ovalnh5[idxm,!istex] <- FALSE
    whm <- which(ovalnh5[idxm,!istex], arr.ind = TRUE)
    whmdf <- cbind(which(idxm)[whm[,"row"]], which(!istex)[whm[,"col"]])
    ovalnh5[whmdf] <- FALSE
  }
  
  ovalnh5
}


## private function .countMultiReadsGenes()
## Counts multi-mapping reads mapping to multiple genes by counting fraction
## counts
.h5_countMultiReadsGenes <- function(ttpar, ovalnh5, maskuniqaln, mt, iste,
                                  istex, tx_idx, readids, alnreadids, ov,
                                  uniqcnt) {
  ovalnh5_multig <- ovalnh5[!maskuniqaln[mt], !istex]
  yesg <- rowSums2(ovalnh5_multig)>0
  
  ## Computing counts for reads with multiple alignments mapping to different
  ## reads and also for multimapping reads with only 1 alignment mapping to
  ## a gene
  ovalnh5_multig <- ovalnh5_multig[yesg,,drop=FALSE]
  
  # Getting the num of different alignments mapping to a gene for each read
  alnreadids_multig <-alnreadids[unique(queryHits(
    ov[!iste[subjectHits(ov)]]))]
  nalnperread <- table(alnreadids_multig) # getting only overlaps from TEs
  readids_multig <- readids[!maskuniqaln[mt]][yesg]
  mt_multig <- match(readids_multig, names(nalnperread))
  
  # Counts provided by unique reads for genes to which multimapping 
  # reads map to
  matmultiguniqc <- t(t(ovalnh5_multig)*uniqcnt[tx_idx][!istex])
  rsum <- rowSums2(matmultiguniqc)
  rsum0 <- rsum == 0
  
  # Reads for which the genes to which the read aligns have > counts
  matmultiguniqc[!rsum0,] <- matmultiguniqc[!rsum0,]/rsum[!rsum0]
  
  # Reads for which the genes to which the read aligns have 0 counts
  matmultiguniqc[rsum0,] <- (ovalnh5_multig[rsum0,] /
                               rowSums(ovalnh5_multig[rsum0,]))
  
  # Adjusting for number of alignments
  matmultiguniqc <- matmultiguniqc/as.numeric(nalnperread[mt_multig])
  
  multigcnt <- rep(0L, length(ttpar@features))
  multigcnt[tx_idx][!istex] <- colSums2(matmultiguniqc)
  multigcnt
}

## private function .h5_countUniqueRead()
## Counts unique reads mapping to TEs and genes (if present)
.h5_countUniqueRead <- function(ttpar, ovalnh5, maskuniqaln, mt, tx_idx, istex) {
  uniqcnt <- rep(0L, length(ttpar@features))
  ovalnh5uniq_g <- ovalnh5[maskuniqaln[mt], !istex, drop=FALSE]
  ovalnh5uniq_te <- ovalnh5[maskuniqaln[mt], istex, drop=FALSE]
  ovmultig <- rowSums2(ovalnh5uniq_g) > 1
  ovmultite <- rowSums2(ovalnh5uniq_te) > 1
  
  # Addressing gene counts: count divided by the number of different genes
  if (any(ovmultig)) {
    # Counting reads overlapping only 1 element
    uniqcnt[tx_idx][!istex] <- colSums2(ovalnh5uniq_g[!ovmultig,])
    # Counting reads overlapping more than 1 element
    mg <- ovalnh5uniq_g[ovmultig,]/rowSums2(ovalnh5uniq_g[ovmultig,])
    uniqcnt[tx_idx][!istex] <- uniqcnt[tx_idx][!istex] + colSums2(mg)
    
  } else {
    uniqcnt[tx_idx][!istex] <- colSums2(ovalnh5uniq_g)
  }
  
  # Addressing TE counts: count divided by the number of different TEs
  # proportionally to the expression level of each TE provided by unique counts
  if (any(ovmultite)) {
    # Counting reads overlapping only 1 element
    uniqcnt[tx_idx][istex] <- colSums2(ovalnh5uniq_te[!ovmultite,])
    # Counting reads overlapping more than 1 element
    mte <- t(t(
      ovalnh5uniq_te[ovmultite,,drop=FALSE])*uniqcnt[tx_idx][istex])
    nocounts <- which(rowSums2(mte) == 0)
    mte[nocounts,,drop=FALSE] <- ovalnh5uniq_te[
      ovmultite,,drop=FALSE][nocounts,]
    mte <- mte/rowSums2(mte)
    uniqcnt[tx_idx][istex] <- uniqcnt[tx_idx][istex] + colSums2(mte)
    
  } else {
    uniqcnt[tx_idx][istex] <- colSums2(ovalnh5uniq_te)
  }
  
  uniqcnt
}


