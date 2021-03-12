#' @importFrom S4Vectors nrow rownames
.checkPhenodata <- function(pdata, nr) {
  if (!is.null(pdata)) {
    if (nrow(pdata) != nr)
      stop("number of rows in 'phenodata' is different than the number of input BAM files in the input parameter object 'x'.")
    if (is.null(rownames(pdata)))
      stop("'phenodata' has no row names.")
  }
}

#' @importFrom S4Vectors DataFrame
.createColumnData <- function(m, pdata) {
  colData <- DataFrame(row.names=gsub(".bam$", "", colnames(m)))
  if (!is.null(pdata))
    colData <- phenodata

  colData
}

## borrowed from GenomicAlignments/R/summarizeOverlaps-methods.R
.checkBamFileListArgs <- function(bfl, singleEnd, fragments) {
  if (singleEnd) {
    if (all(isTRUE(asMates(bfl))))
      stop("cannot specify both 'singleEnd=TRUE' and 'asMates=TRUE'")
    if (fragments)
      stop("when 'fragments=TRUE', 'singleEnd' should be FALSE")
  } else
    asMates(bfl) <- TRUE
}

## borrowed from GenomicAlignments/R/summarizeOverlaps-methods.R
.getReadFunction <- function(singleEnd, fragments) {
  if (singleEnd) {
    FUN <- readGAlignments
  } else {
    if (fragments)
      FUN <- readGAlignmentsList
    else
      FUN <- readGAlignmentPairs
  }
  
  FUN
}

#' @importFrom S4Vectors nLnode nRnode isSorted from to Hits
## appends the second Hits object to the end of the first one
## assuming they have identical right nodes
.appendHits <- function(hits1, hits2) {
  stopifnot(nRnode(hits1) == nRnode(hits2))
  stopifnot(isSorted(from(hits1)) == isSorted(from(hits2)))
  hits <- c(Hits(from=from(hits1), to=to(hits1),
                 nLnode=nLnode(hits1)+nLnode(hits2),
                 nRnode=nRnode(hits1), sort.by.query=isSorted(from(hits1))),
            Hits(from=from(hits2)+nLnode(hits1), to=to(hits2),
                 nLnode=nLnode(hits1)+nLnode(hits2),
                 nRnode=nRnode(hits2), sort.by.query=isSorted(from(hits2))))
  hits
}

