## private function .checkPhenoData()

#' @importFrom S4Vectors nrow rownames
.checkPhenodata <- function(pdata, nr) {
  if (!is.null(pdata)) {
    if (nrow(pdata) != nr)
      stop("number of rows in 'phenodata' is different than the number of input BAM files in the input parameter object 'x'.")
    if (is.null(rownames(pdata)))
      stop("'phenodata' has no row names.")
  }
}

## private function .createColumnData()

#' @importFrom S4Vectors DataFrame
.createColumnData <- function(m, pdata) {
  colData <- DataFrame(row.names=gsub(".bam$", "", colnames(m)))
  if (!is.null(pdata))
    colData <- pdata

  colData
}

## private function .checkBamFileListArgs()
## adapted from GenomicAlignments/R/summarizeOverlaps-methods.R

#' @importFrom Rsamtools BamFileList asMates
.checkBamFileListArgs <- function(bfl, singleEnd, fragments) {
  if (missing(bfl) || !class(bfl) %in% c("character", "BamFileList"))
    stop("argument 'bfl' should be either a string character vector of BAM file names or a 'BamFileList' object")

  if (is.character(bfl)) {
    mask <- sapply(bfl, file.exists)
    if (any(!mask))
      stop(sprintf("The following input BAM files cannot be found:\n%s",
                   paste(paste("  ", bfl), collapse="\n")))
  }

  if (!is(bfl, "BamFileList"))
    bfl <- BamFileList(bfl, asMates=!singleEnd)

  if (singleEnd) {
    if (all(isTRUE(asMates(bfl))))
      stop("cannot specify both 'singleEnd=TRUE' and 'asMates=TRUE'")
    if (fragments)
      stop("when 'fragments=TRUE', 'singleEnd' should be FALSE")
  } else
    asMates(bfl) <- TRUE

  bfl
}

## private function .processFeatures()
## builds a single 'GRanges' object from input TE and gene features.
## parameters: teFeatures - a 'GRanges' or 'GRangesList' object with
##                          TE annotations
##             teFeaturesobjname - the name of 'teFeatures'
##             geneFeatures - a 'GRanges' or 'GRangesList' object with
##                            gene annotations
##             geneFeaturesobjname - the name of 'geneFeatures'
##             aggregateby - names of metadata columns in 'teFeatures'
##                           to be used later for aggregating estimated
##                           counts.

#' @importFrom S4Vectors mcols Rle
.processFeatures <- function(teFeatures, teFeaturesobjname, geneFeatures,
                             geneFeaturesobjname, aggregateby) {

  if (!exists(teFeaturesobjname))
    stop(sprintf("input TE features object '%s' is not defined.",
                 teFeaturesobjname))

  if (!is(teFeatures, "GRanges") && !is(teFeatures, "GRangesList"))
    stop(sprintf("TE features object '%s' should be either a 'GRanges' or a 'GRangesList' object.",
                 teFeaturesobjname))

  if (length(aggregateby) > 0)
    if (any(!aggregateby %in% colnames(mcols(teFeatures))))
        stop(sprintf("%s not in metadata columns of the TE features object.",
             paste(aggregateby[!aggregateby %in% colnames(mcols(teFeatures))])))

  if (is.null(names(teFeatures)) && length(aggregateby) == 0)
    stop(sprintf("the TE features object '%s' has no names and no aggregation metadata columns have been specified.",
                 teFeaturesobjname))

  features <- teFeatures
  if (is(teFeatures, "GRangesList"))
    features <- unlist(teFeatures)

  if (!is.na(geneFeatures)) {
    geneFeaturesobjname <- deparse(substitute(geneFeatures))
    if (!is(geneFeatures, "GRanges") && !is(geneFeatures, "GRangesList"))
      stop(sprintf("gene features object '%s' should be either a 'GRanges' or a 'GRangesList' object.",
                   geneFeaturesobjname))
    if (any(names(geneFeatures) %in% names(teFeatures)))
      stop("gene features have some common identifiers with the TE features.")

    if (length(geneFeatures) == 0)
      stop(sprintf("gene features object '%s' is empty.", geneFeaturesobjname))

    if (is(geneFeatures, "GRangesList"))
      geneFeatures <- unlist(geneFeatures)

    features <- c(teFeatures, geneFeatures)
    temask <- Rle(rep(FALSE, length(features) + length(geneFeatures)))
    temask[1:length(features)] <- TRUE
    features$isTE <- temask
  }

  features
}

## private function .consolidateFeatures()
## builds a 'GRanges' or 'GRangesList' object
## grouping TE features, if necessary, and
## adding gene features, if available.
## parameters: x - TEtranscriptsParam object
##             fnames - feature names vector to which
##                      consolidated features should match

#' @importFrom methods is
#' @importFrom S4Vectors split
#' @importFrom GenomicRanges GRangesList
.consolidateFeatures <- function(x, fnames) {
  teFeatures <- x@annotations
  if (!is.null(x@annotations$isTE) && any(x@annotations$isTE)) {
    teFeatures <- x@annotations[x@annotations$isTE]
  }

  if (length(x@aggregateby) > 0) {
    f <- .factoraggregateby(teFeatures, x@aggregateby)
    teFeatures <- split(teFeatures, f)
  }

  features <- teFeatures
  if (!is.null(x@annotations$isTE) && any(!x@annotations$isTE)) {
    geneFeatures <- x@annotations[!x@annotations$isTE]
    if (is(features, "GRangesList")) ## otherwise is a GRanges object
      geneFeatures <- split(geneFeatures, names(geneFeatures))
    features <- c(features, geneFeatures)
  }

  stopifnot(length(features) == length(fnames)) ## QC
  features <- features[match(fnames, names(features))]

  features
}

## private function .factoraggregateby()
## builds a factor with as many values as the
## length of the input annotations in 'ann', where
## every value is made by pasting the columns in
## 'aggby' separated by ':'.
## parameters: ann - GRanges object with annotations
##             aggby - names of metadata columns in 'ann'
##                     to be pasted together

#' @importFrom GenomicRanges mcols
f <- .factoraggregateby <- function(ann, aggby) {
  stopifnot(all(aggby %in% colnames(mcols(ann)))) ## QC
  spfstr <- paste(rep("%s", length(aggby)), collapse=":")
  f <- do.call("sprintf", c(spfstr, as.list(mcols(ann)[, aggby])))
  f
}

## private function .getReadFunction()
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

## private function .appendHits()
## appends the second Hits object to the end of the first one
## assuming they have identical right nodes

#' @importFrom S4Vectors nLnode nRnode isSorted from to Hits
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
