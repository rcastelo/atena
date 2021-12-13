## pretty-print names (private)
.pprintnames <- function(x) {
    y <- x
    if (length(x) > 2)
        y <- c(y[1], "...", y[length(y)])
    y <- paste(y, collapse=", ")
    y
}

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

#' @importFrom Rsamtools BamFileList asMates asMates<-
.checkBamFileListArgs <- function(bfl, singleEnd, fragments) {
    if (missing(bfl) || !class(bfl) %in% c("character", "BamFileList"))
        stop("argument 'bfl' should be either a string character vector of BAM file names or a 'BamFileList' object")
    
    if (is.character(bfl)) {
        mask <- vapply(bfl, FUN = file.exists, FUN.VALUE = logical(1))
        if (any(!mask))
            stop(sprintf("The following input BAM files cannot be found:\n%s",
                        paste(paste("  ", bfl[!mask]), collapse="\n")))
    }
    
    if (!is(bfl, "BamFileList"))
        bfl <- BamFileList(bfl, asMates=!singleEnd)
    
    if (singleEnd) {
        if (all(isTRUE(asMates(bfl))))
            stop("cannot specify both 'singleEnd=TRUE' and 'asMates=TRUE'")
        # if (fragments)
        #     stop("when 'fragments=TRUE', 'singleEnd' should be FALSE")
    } else
        asMates(bfl) <- TRUE
    
    bfl
}

## private function .checkBamReadMapper()
## extracts the name of the read mapper software from one or more BAM files
## parameters: bamfiles - BAM file names

#' @importFrom Rsamtools scanBamHeader
.checkBamReadMapper <- function(bamfiles) {
    if (missing(bamfiles) || !"character" %in% class(bamfiles))
        stop("argument 'bamfiles' should be a string character vector of BAM file names")
    
    mask <- vapply(bamfiles, FUN = file.exists, FUN.VALUE = logical(1L))
    if (any(!mask))
        stop(sprintf("The following input BAM files cannot be found:\n%s",
                    paste(paste("  ", bamfiles[!mask]), collapse="\n")))
    
    hdr <- scanBamHeader(bamfiles)
    readaligner <- vapply(hdr, FUN = function(x) {
                            ra <- NA_character_
                            if (!is.null(x$text[["@PG"]])) {
                                pgstr <- x$text[["@PG"]]
                                mt <- gregexpr("^PN:", pgstr)
                                wh <- which(vapply(mt, FUN = function(x) x!=-1,
                                            FUN.VALUE = logical(1L)))
                                ra <- substr(pgstr[[wh]],
                                            attr(mt[[wh]], "match.length") + 1,
                                            100000L)
                            }
                            tolower(ra)
                    }, FUN.VALUE = character(1))
    readaligner <- readaligner[!duplicated(readaligner)]
    readaligner <- as.vector(readaligner[!is.na(readaligner)])
    if (length(readaligner) == 0)
        warning("no read aligner software information in BAM files.")
    if (any(readaligner[1] != readaligner))
        warning(sprintf("different read aligner information in BAM files. Assuming %s",
                        readaligner[1]))
    
    readaligner[1]
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

#' @importFrom S4Vectors mcols Rle DataFrame
#' @importFrom GenomeInfoDb seqlevels<- seqlevels
.processFeatures <- function(teFeatures, teFeaturesobjname, geneFeatures,
                                geneFeaturesobjname, aggregateby,
                                aggregateexons) {
    
    if (missing(teFeatures)) 
        stop("missing 'teFeatures' argument.")
    
    # if (!exists(teFeaturesobjname))
    #     stop(sprintf("input TE features object '%s' is not defined.",
    #                 teFeaturesobjname))
    
    if (!is(teFeatures, "GRanges") && !is(teFeatures, "GRangesList"))
        stop(sprintf("TE features object '%s' should be either a 'GRanges' or a 'GRangesList' object.",
                    teFeaturesobjname))
    
    if (length(aggregateby) > 0)
        if (any(!aggregateby %in% colnames(mcols(teFeatures))))
            stop(sprintf("%s not in metadata columns of the TE features object.",
                aggregateby[!aggregateby %in% colnames(mcols(teFeatures))]))
    
    if (is.null(names(teFeatures)) && length(aggregateby) == 0)
        stop(sprintf("the TE features object '%s' has no names and no aggregation metadata columns have been specified.",
                    teFeaturesobjname))
    
    features <- teFeatures
    if (is(teFeatures, "GRangesList")) 
        features <- unlist(teFeatures)
    
    if (!all(is.na(geneFeatures))) {
        if (is(geneFeatures, "GRangesList")) 
            geneFeatures <- unlist(geneFeatures)
        
        features <- .joinTEsGenes(teFeatures, geneFeatures)
    } else {
        features$isTE <- rep(TRUE, length(features))
    }
    
    iste <- as.vector(features$isTE)
    if (!all(is.na(geneFeatures))) {
        if (aggregateexons & !all(iste) & !is.null(mcols(geneFeatures)$type)) {
            iste <- aggregate(iste, by = list(names(features)), unique)
            features <- .groupGeneExons(features)
            mtname <- match(names(features), iste$Group.1)
            iste <- iste[mtname,"x"]
        }
    }
    
    attr(features, "isTE") <- DataFrame("isTE" = iste)
    features
}


## private function .groupGeneExons()
## groups exons from the same gene creating a 'GRangesList' object
.groupGeneExons <- function(features) {
    if (!any(mcols(features)$type == "exon")) {
        stop(".groupGeneExons: no elements with value 'exon' in 'type' column of the metadata of the 'GRanges' or 'GRangesList' object with gene annotations.")
    }
    features <- features[features$isTE | features$type == "exon"]
    featuressplit <- split(x = features, f = names(features))
    featuressplit
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
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges IRanges
.consolidateFeatures <- function(x, fnames) {
    
    iste <- as.vector(attributes(x@features)$isTE[,1])
    teFeatures <- x@features
    if (!is.null(iste) && any(iste)) {
        teFeatures <- x@features[iste]
    }
    
    if (length(x@aggregateby) > 0) {
        f <- .factoraggregateby(teFeatures, x@aggregateby)
        if (is(teFeatures, "GRangesList"))
            teFeatures <- unlist(teFeatures)
        teFeatures <- split(teFeatures, f)
    }
    
    features <- teFeatures
    if (!is.null(iste) && any(!iste)) {
        geneFeatures <- x@features[!iste]
        # if (is(features, "GRangesList")) ## otherwise is a GRanges object
        #   geneFeatures <- split(geneFeatures, names(geneFeatures))
        features <- c(features, geneFeatures)
    }
    
    stopifnot(length(features) == length(fnames)) ## QC
    features <- features[match(fnames, names(features))]
    
    if (is(x, "TelescopeParam")) {
        nofeat_gr <- GRanges(seqnames = "chrNofeature", 
                             ranges = IRanges(start = 1, end = 1))
        names(nofeat_gr) <- "no_feature"
        seqlev <- unique(c(seqlevels(features),seqlevels(nofeat_gr)))
        seqlevels(features) <- seqlev
        seqlevels(nofeat_gr) <- seqlev
        features <- c(features, nofeat_gr)
    }
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
.factoraggregateby <- function(ann, aggby) {
    if (is(ann,"GRangesList")) {
        ann <- unlist(ann)
    }
    stopifnot(all(aggby %in% colnames(mcols(ann)))) ## QC
    if (length(aggby) == 1) {
        f <- mcols(ann)[, aggby]
    } else {
        spfstr <- paste(rep("%s", length(aggby)), collapse=":")
        f <- do.call("sprintf", c(spfstr, as.list(mcols(ann)[, aggby])))
    }
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


#' @importFrom GenomeInfoDb seqlevels<- seqlevels
.joinTEsGenes <- function(teFeatures, geneFeatures) {
    
    geneFeaturesobjname <- deparse(substitute(geneFeatures))
    if (!is(geneFeatures, "GRanges") && !is(geneFeatures, "GRangesList"))
        stop(sprintf("gene features object '%s' should be either a 'GRanges' or a 'GRangesList' object.",
                    geneFeaturesobjname))
    if (any(names(geneFeatures) %in% names(teFeatures)))
        stop("gene features have some common identifiers with the TE features.")
    
    if (length(geneFeatures) == 0)
        stop(sprintf("gene features object '%s' is empty.", geneFeaturesobjname))
    
    slev <- unique(c(seqlevels(teFeatures), seqlevels(geneFeatures)))
    seqlevels(teFeatures) <- slev
    seqlevels(geneFeatures) <- slev
    features <- c(teFeatures, geneFeatures)
    temask <- Rle(rep(FALSE, length(teFeatures) + length(geneFeatures)))
    temask[seq_along(teFeatures)] <- TRUE
    features$isTE <- temask
    features
}

#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix Matrix
.buildOvValuesMatrix <- function(x, ov, values, aridx, ridx, fidx) {
    stopifnot(class(values) %in% c("logical", "integer", "numeric")) ## QC
    ovmat <- Matrix(do.call(class(values), list(1)),
                    nrow=length(ridx), ncol=length(fidx))
    mt1 <- match(aridx[queryHits(ov)], ridx)
    mt2 <- match(subjectHits(ov), fidx)
    mtov <- cbind(mt1, mt2)
    if (is(x, "TelescopeParam")) {
        # mtov <- cbind(mt1, mt2)
        mtalign <- match(paste(mtov[,1],mtov[,2],sep = ":"),
                         unique(paste(mtov[,1],mtov[,2], sep = ":")))
        s <- split(x = values[queryHits(ov)], f = mtalign)
        smax <- unlist(lapply(s, max), use.names = FALSE)
        values <- smax[mtalign]
    } else {
        values <- values[queryHits(ov)]
    }
    
    ovmat[mtov] <- values
    ovmat
}
