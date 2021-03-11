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
