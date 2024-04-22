#' Get RepeatMasker UCSC annotations
#'
#' The \code{annotaTEs()} function fetches RepeatMasker UCSC transposable
#' element (TE) annotations using 
#' \link[AnnotationHub]{AnnotationHub} and parses them.
#' 
#' @param genome The genome version of the desired RepeatMasker annotations
#'        (e.g. "hg38").
#'
#' @param parsefun A \code{function} to parse the annotations:
#' \itemize{
#'   \item Function \code{rmskidentity} returns RepeatMasker annotations as
#'         present in \link[AnnotationHub]{AnnotationHub}, without
#'         processing them.
#'   \item Function \code{rmskbasicparser} parses annotations by removing
#'         low complexity regions, simple repeats, satellites, rRNA, scRNA,
#'         snRNA, srpRNA and tRNA. Also removes TEs
#'         with a strand different than "+" or "-". Modifies "repFamily" and
#'         "repClass" columns when a "?" is present or when they are defined
#'         as "Unknown" or "Other". Finally, assigns a unique id to each TE
#'         instance by adding the suffix "_dup" plus a number at the end of
#'         the "repName".
#'   \item Function \code{rmskatenaparser} parses RepeatMasker annotations 
#'         reconstructing fragmented TEs by assembling together fragments from
#'         the same TE that are close enough. For LTR class TEs, it tries to
#'         reconstruct full-length and partial TEs following the LTR - internal
#'         region - LTR structure. Input is a \code{GRanges} object and output
#'         is a \code{GRangesList} object.
#'   \item Function \code{OneCodeToFindThemAll} parses annotations following
#'         the 'One code to find them all' method by 
#'         \href{https://doi.org/10.1186/1759-8753-5-13}{(Bailly-Bechet et al. 2014)}. 
#'         Input is a \code{GRanges} object and output is a \code{GRangesList} 
#'         object.
#'   \item User-defined function. Input and output should be \code{GRanges}
#'         objects.
#' }
#' 
#' @param AHid AnnotationHub unique identifier, of the form AH12345, of an
#'             object with TE annotations. This is an optional argument to
#'             specify a concrete AnnotationHub resource, for instance
#'             when more there is more than one RepeatMasker annotation
#'             available for a specific genome version. If \code{AHid} is
#'             not specified, the latest RepeatMasker annotation is be used.
#' 
#' @param ... Arguments passed to \code{parsefun}.
#'        
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object with
#'         transposable element annotations.
#'
#' @details
#' Given a specific genome version, the \code{annotaTEs()} function fetches
#' RepeatMasker annotations from UCSC Genome Browser using the 
#' \link[AnnotationHub]{AnnotationHub} package. Since RepeatMasker not only
#' provides TE annotations but also low complexity DNA sequences and other
#' types of repeats, a specific \code{parsefun} can be set to parse these
#' annotations (e.g. \code{rmskbasicparser} or a user-defined function). If no
#' parsing is required, \code{parsefun} can be set to \code{rmskidentity}.
#'
#' @seealso
#' \code{\link[AnnotationHub]{AnnotationHub}}
#' 
#'
#' @examples
#' rmskid <- annotaTEs(genome="hg19", parsefun=rmskidentity)
#' 
#' 
#' @aliases annotaTEs
#' @rdname annotaTEs
#' @name annotaTEs
#' @importFrom AnnotationHub AnnotationHub query
#' @export
annotaTEs <- function(genome="hg38", parsefun=rmskidentity, AHid = NULL, ...) {
    suppressMessages(ah <- AnnotationHub())
    if (!is.null(AHid))
        qah <- ah[AHid]
    else
        suppressMessages(qah <- query(ah, c(genome, "RepeatMasker", "UCSC")))
    
    if (length(qah) == 0) {
        stop(sprintf("UCSC RepeatMasker tracks for genome %s not found", 
                     genome))
    }
    else if (length(qah) > 1) {
        message(sprintf("more than one UCSC RepeatMasker track for genome %s found, using the latest one", genome))
        mt <- gregexpr(pattern=" \\([A-Za-z0-9]+\\) ", qah$title)
        qahdates <- substr(qah$title, unlist(mt)+2,
                             unlist(mt)+sapply(mt, attr, "match.length")-3)
        qahdates <- as.Date(paste0("1", qahdates), "%d%b%Y")
        ## in case > 1 RM annotations available, pick the most recent one
        suppressMessages(qah <- qah[which.max(qahdates)])
    }
    
    id <- names(qah)
    gr <- ah[[id]]
    parsefun(gr, ...)
}

#' Parser of RepeatMasker annotations
#' @param gr A \link[GenomicRanges:GRanges-class]{GRanges} object with
#'           RepeatMasker annotations from \link[AnnotationHub]{AnnotationHub}
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object.
#'         
#' @details 
#' Parses annotations by removing low complexity regions, simple repeats,
#' satellites, rRNA, scRNA, snRNA, srpRNA and tRNA. Also removes TEs with a
#' strand different than "+" or "-". Modifies "repFamily" and
#' "repClass" columns when a "?" is present or when they are defined
#' as "Unknown" or "Other". Finally, assigns a unique id to each TE instance
#' by adding the suffix "_dup" plus a number at the end of the "repName".
#'
#' @examples
#' rmskba <- annotaTEs(genome="dm6", parsefun=rmskbasicparser)
#' 
#' @aliases rmskbasicparser
#' @rdname rmskbasicparser
#' @name rmskbasicparser
#' @importFrom GenomicRanges strand
#' @export
rmskbasicparser <- function(gr) {
    to_discard <- c("Low_complexity", "Simple_repeat", "Satellite", "rRNA",
                    "scRNA", "snRNA", "srpRNA", "tRNA")
    gr <- gr[!(gr$repClass %in% to_discard),]
    gr <- gr[strand(gr) == "+" | strand(gr) == "-",]
    # removing '?' from name, class or family
    gr$repName <- gsub(pattern="?", x=gr$repName, replacement="", fixed =TRUE)
    gr$repClass <- gsub(pattern="?", x=gr$repClass, replacement="", fixed=TRUE)
    gr$repFamily <- gsub(pattern="?", x=gr$repFamily, replacement="", 
                         fixed =TRUE)
    gr <- unlist(split(gr, f=gr$repName))
    dupl_te <- duplicated(gr$repName) | duplicated(gr$repName, fromLast = TRUE)
    ndup <- lengths(split(gr, f=gr$repName))
    dup <- unlist(lapply(ndup, FUN=function(x) seq_len(x)))
    repName <- paste(gr$repName, "_dup", dup, sep="")
    repName[!dupl_te] <- gr$repName[!dupl_te]
    names(gr) <- repName
    gr <- sort(gr)

    i <- gr$repFamily %in% c("Unknown", "Other")
    if (any(i)) {
      gr$repFamily[i] <- gr$repName[i]
    }
    gr
}

#' Identity function for parsefun
#' @param gr A \link[GenomicRanges:GRanges-class]{GRanges} object.
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object.
#'
#' @details 
#' Identity function: returns the \link[GenomicRanges:GRanges-class]{GRanges}
#' object without any modification.
#' 
#' @examples
#' rmskid <- annotaTEs(genome="dm6", parsefun=rmskidentity)
#'         
#' @aliases rmskidentity
#' @rdname rmskidentity
#' @name rmskidentity
#' @export
rmskidentity <- function(gr) {
  gr
}

#' OneCodeToFindThemAll parser of RepeatMasker annotations
#' @param gr A \link[GenomicRanges:GRanges-class]{GRanges} object with
#' RepeatMasker annotations from \link[AnnotationHub]{AnnotationHub}
#' 
#' @param dictionary (Default NULL) When NULL, a dictionary is built based 
#' on names of repeats. If not, a data.frame with equivalences LTR - internal
#' regions created by the user, where first column should be the name of the
#' internal region and the second column should be the LTR(s). When more than
#' one LTR, these should be separated by ":".
#' 
#' @param fuzzy (Default FALSE) A logical; if TRUE, the search for equivalences
#' between internal parts and LTRs to reconstruct LTR class transposable 
#' elements is less stringent, allowing more matches between corresponding
#' subparts. This option can increase the proportion of false positives
#' (incorrectly reconstructed LTR class TEs).
#' 
#' @param strict (Default FALSE) A logical; if TRUE, the 80-80 rule is applied,
#' i.e. only copies with more than 80% identity to the reference
#' and more than 80 bp long are reported.
#' 
#' @param insert (Default -1) An integer. When \code{insert} < 0, two fragments
#' are assembled if the distance separating their furthest extremities is less 
#' than twice the reference length of the element. When \code{insert} > 0,
#' fragments are assembled if the distance between their closest extremities
#' is equal or less than \code{insert}. When \code{insert} = 0, two fragments
#' are assembled if they are in contact next to each other.
#' 
#' @param BPPARAM See \code{?\link[BiocParallel:bplapply]{bplapply}} in the 
#' BiocParallel package. Can be used to run function in parallel.
#'
#' @return A \link[GenomicRanges:GRangesList-class]{GRangesList} object.
#'         
#' @details 
#' Implementation of One code to find them all 
#' \href{https://doi.org/10.1186/1759-8753-5-13}{(Bailly-Bechet et al. 2014)}.
#' Parses RepeatMasker annotations from UCSC by assembling together fragments
#' from the same transposable elemenet (TE) that are close enough (determined
#' by the \code{insert} parameter). For TEs from the LTR class, the parser
#' tries to reconstruct full-length, when possible, or partial TEs following
#' the LTR - internal region - LTR structure. Equivalences between internal
#' regions and flanking LTRs can be set by the user with the \code{dictionary}
#' parameter or can be obtained by the parser. In this last case, the 
#' \code{fuzzy} parameter determines the level of stringency when searching 
#' for LTR - internal region equivalences.
#'
#' @examples
#' \dontrun{
#' rmskoc <- annotaTEs(genome="dm6", parsefun=OneCodeToFindThemAll,
#'                     fuzzy=FALSE, strict=FALSE)
#' }
#' 
#' @references
#' Bailly-Bechet et al. "One code to find them all": a perl tool to 
#' conveniently parse RepeatMasker output files. Mobile DNA. 2014;5(1):1-15.
#' DOI: \url{https://doi.org/10.1186/1759-8753-5-13}
#' 
#' @aliases OneCodeToFindThemAll
#' @rdname OneCodeToFindThemAll
#' @name OneCodeToFindThemAll
#' @importFrom GenomicRanges strand width mcols "mcols<-"
#' @importFrom BiocParallel bplapply
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors runValue
#' @importFrom IRanges mean
#' @export
OneCodeToFindThemAll <- function(gr, dictionary=NULL, fuzzy = FALSE,
                                 strict = FALSE, insert = -1, 
                                 BPPARAM = SerialParam(progressbar=TRUE)) {

  if (!is.integer(insert) & !is.numeric(insert))
    stop("'insert' must be an integer value")
  
  if (!is(gr, "GRanges"))
    stop("'gr' should be a GRanges object with RepeatMasker annotations")
  
  if (length(gr) == 0)
    stop("'gr' is empty")
  
  if (!is.null(dictionary)) {
    if (!file.exists(dictionary))
        stop("'dictionary' should be NULL or specify the name of a dictionary",
        " file, see ?OneCodeToFindThemAll")
  }
  
  if (is.null(dictionary)) {
    inout <- .builDictionary(gr, fuzzy)
  } else {
    inside <- dictionary[, 2]
    names(inside) <- dictionary[, 1]
    outsidesp <- strsplit(dictionary[, 2], split = ":")
    outside <- rep(dictionary[, 1], lengths(outsidesp))
    names(outside) <- unlist(outsidesp)
    inout <- list(inside = inside, outside = outside)
  }
  
  gr <- .filterNonTEs(gr)
  gr$repStart <- abs(gr$repStart)
  gr$repLeft <- abs(gr$repLeft)
  cons_length <- .getConsLength(gr)
  
  # Parsing by chromosome
  annrec <- bplapply(runValue(seqnames(gr)), .OCTFTA_parser, gr=gr, 
                     cons_length=cons_length, outside=inout$outside, 
                     inside=inout$inside, insert=insert, BPPARAM=BPPARAM)
  
  # 'annrec' is a list of GRangesLists, here we concatenate GRangesList to a
  # single GRangesList
  annrec <- do.call("c", annrec)
  
  # -- Filtering based on 'strict' option
  if (strict) {
    maxdiv <- mean(relist(mcols(unlist(annrec))$milliDiv, annrec)) < 200
    tokeep <- sum(width(annrec)) >= 80 & maxdiv
    annrec <- annrec[tokeep]
  }
  
  # -- Computation relative length --
  # for full-length and partial ERVs, and ERVs with only internal region, the
  # consensus length is considered that of the full-length TE: LTR + int + LTR
  nTEs <- suppressWarnings(do.call("rbind", strsplit(names(annrec), ".", 
                                                     fixed=TRUE)))[,1]
  mcols(annrec)$RelLength <- .getRelLength(nTEs, cons_length, inout, annrec)
  
  # -- Simplifying TE name: subfamily name + number of TE separated by "." --
  names(annrec) <- paste(nTEs, seq_along(annrec), sep = ".")
  
  # -- Including TE class name in GRangesList 
  mcols(annrec)$Class <- unlist(unique(relist(unlist(annrec)$repClass, annrec)))
  
  annrec
}

#' Getter functions of TE classes from parsed RepeatMasker annotations.
#'
#' @param annot A [`GRanges`] or [`GRangesList`] object obtained with the
#'              function `annotaTES()`, using either [`OneCodeToFindThemAll`]
#'              or [`rmskatenaparser`] as RepeatMasker parser functions.
#'              Alternatively, if `annot` is a [`QuantifyParam`] or a
#'              [`SummarizedExperiment`] object produced by the `qtex()`
#'              function, this function will attempt to extract the
#'              corresponding annotations from inside those objects.
#' 
#' @param relLength (Default 0.9) Numeric value that can take values between 0
#'                  to 1. Sets the minimum relative length required for
#'                  features. Elements with a lower relative length than
#'                  \code{relLength} will be filtered. The relative length
#'                  used is the one obtained by \code{OneCodeToFindThemAll()}
#'                  or \code{rmskatenaparser()}.
#'                  (length of the reconstructed TE / length of the reference).
#' 
#' @param fullLength (Default TRUE) Logical value on whether reconstructed
#'                   full-length LTR TEs (elements with structure 
#'                   LTR - internal region - LTR) should be selected.
#' 
#' @param partial (Default FALSE) Logical value on whether partially
#'                reconstructed LTR TEs should be selected (structure LTR -
#'                internal region or internal region - LTR).
#'                               
#' @param soloLTR (Default FALSE) Logical value on whether solo LTRs should be
#'                selected. Note that only fragments unambiguously identified
#'                as LTRs thanks to the identification of their equivalent
#'                internal region are considered as LTRs.
#' 
#' @param otherLTR (Default FALSE) Logical value on whether other TEs from the
#'                 LTR class, not included in any of the previous three
#'                 categories, should be selected. These include TEs from LTR
#'                 class that cannot be unambiguously identified as LTR o
#'                 internal region, and thus cannot be reconstructed into
#'                 partial or full-length elements; as well as solo internal
#'                 regions. 
#'
#' @param returnMask (Default FALSE) Logical value indicating whether a subset
#'                   of the input annotations should be returned (default) or
#'                   a logical mask of the same length as the input annotations
#'                   where \code{TRUE} values indicate what annotations belong
#'                   to the TE class we want to obtain with the getter function.
#'
#' @return A \link[GenomicRanges:GRangesList-class]{GRangesList} object with
#'         annotations from class corresponding to the getter function (LTRs,
#'         LINEs, SINEs or DNA transposons).
#'         
#' @details 
#' Retrieves annotations from the TE class corresponding to the getter function,
#' using RepeatMasker annotations after parsing them with the
#' \code{OneCodeToFindThemAll()} or \code{rmskatenaparser()} function. The
#' \code{relLength} parameter can be used to filter out elements with a lower
#' relative length. Further parameters can be used to fine-tune the type of
#' elements to be reported.
#'
#' @examples
#' rmskat <- annotaTEs(genome="dm6", parsefun=rmskatenaparser,
#'                     strict=FALSE)
#'
#' rmskat_ltr <- getLTRs(rmskat, relLength=0.95, fullLength=TRUE,
#'                       partial=TRUE)
#'  
#' @aliases getLTRs
#' @rdname annotaTEsGetters
#' @name annotateTEsGetters
#' @importFrom GenomicRanges mcols
#' @export
getLTRs <- function(annot, relLength=0.9, fullLength=TRUE, partial=FALSE,
                    soloLTR=FALSE, otherLTR=FALSE, returnMask=FALSE) {
    
    annot <- .fetchAnnotationsFromObject(annot)

    cnames <- c("Status", "RelLength", "Class")
    mask <- cnames %in% colnames(mcols(annot))
    if (any(!mask))
        stop(sprintf("column(s) %s missing from input annotations.",
                     paste(cnames[!mask], collapse=", ")))
    if (!any(c(fullLength, partial, soloLTR, otherLTR)))
        stop("at least one of these arguments should be TRUE: 'fullLength',",
             " 'partial', 'soloLTR', 'otherLTR'")

    mask <- rep(TRUE, length(annot))
    if ("isTE" %in% colnames(mcols(annot))) {
        mask <- mcols(annot)$isTE
        annot <- annot[mcols(annot)$isTE]
    }
  
    keep <- mcols(annot)$Class == "LTR" & 
            mcols(annot)$RelLength > relLength
    LTRtypes <- c("full-lengthLTR"[fullLength], "partialLTR_up"[partial], 
                  "partialLTR_down"[partial], "LTR"[soloLTR],
                  "int"[otherLTR], "noLTR"[otherLTR])
    keep2 <- mcols(annot)$Status %in% LTRtypes

    res <- mask
    if ("isTE" %in% colnames(mcols(annot))) {
        res[mask] <- keep & keep2
        mcols(annot)$isTE <- NULL
    } else
        res <- keep & keep2

    if (!returnMask)
        res <- annot[keep & keep2]

    res
}

#' @examples
#' rmskat_line <- getLINEs(rmskat, relLength=0.95)
#'  
#' @aliases getLINEs
#' @rdname annotaTEsGetters
#' @name annotateTEsGetters
#' @importFrom GenomicRanges mcols
#' @export
getLINEs <- function(annot, relLength=0.9, returnMask=FALSE) {

    annot <- .fetchAnnotationsFromObject(annot)

    cnames <- c("RelLength", "Class")
    mask <- cnames %in% colnames(mcols(annot))
    if (any(!mask))
        stop(sprintf("column(s) %s missing from input annotations.",
                     paste(cnames[!mask], collapse=", ")))

    mask <- rep(TRUE, length(annot))
    if ("isTE" %in% colnames(mcols(annot))) {
        mask <- mcols(annot)$isTE
        annot <- annot[mcols(annot)$isTE]
    }
  
    keep <- mcols(annot)$Class == "LINE" &
            mcols(annot)$RelLength > relLength

    res <- mask
    if ("isTE" %in% colnames(mcols(annot))) {
        res[mask] <- keep
        mcols(annot)$isTE <- NULL
    } else
        res <- keep

    if (!returnMask)
        res <- annot[keep]

    res
}

#' @examples
#' rmskat_sine <- getSINEs(rmskat, relLength=0.95)
#'  
#' @aliases getSINEs
#' @rdname annotaTEsGetters
#' @name annotateTEsGetters
#' @importFrom GenomicRanges mcols
#' @export
getSINEs <- function(annot, relLength=0.9, returnMask=FALSE) {

    annot <- .fetchAnnotationsFromObject(annot)

    cnames <- c("RelLength", "Class")
    mask <- cnames %in% colnames(mcols(annot))
    if (any(!mask))
        stop(sprintf("column(s) %s missing from input annotations.",
                     paste(cnames[!mask], collapse=", ")))

    mask <- rep(TRUE, length(annot))
    if ("isTE" %in% colnames(mcols(annot))) {
        mask <- mcols(annot)$isTE
        annot <- annot[mcols(annot)$isTE]
    }
  
    keep <- mcols(annot)$Class == "SINE" &
            mcols(annot)$RelLength > relLength

    res <- mask
    if ("isTE" %in% colnames(mcols(annot))) {
        res[mask] <- keep
        mcols(annot)$isTE <- NULL
    } else
        res <- keep

    if (!returnMask)
        res <- annot[keep]

    res
}

#' @examples
#' rmskat_DNAtrans <- getDNAtransposons(rmskat, relLength=0.95)
#'  
#' @aliases getDNAtransposons
#' @rdname annotaTEsGetters
#' @name annotateTEsGetters
#' @importFrom GenomicRanges mcols
#' @export
getDNAtransposons <- function(annot, relLength=0.9, returnMask=FALSE) {

    annot <- .fetchAnnotationsFromObject(annot)

    cnames <- c("RelLength", "Class")
    mask <- cnames %in% colnames(mcols(annot))
    if (any(!mask))
        stop(sprintf("column(s) %s missing from input annotations.",
                     paste(cnames[!mask], collapse=", ")))

    mask <- rep(TRUE, length(annot))
    if ("isTE" %in% colnames(mcols(annot))) {
        mask <- mcols(annot)$isTE
        annot <- annot[mcols(annot)$isTE]
    }
  
    keep <- mcols(annot)$Class == "DNA" & 
            mcols(annot)$RelLength > relLength

    res <- mask
    if ("isTE" %in% colnames(mcols(annot))) {
        res[mask] <- keep
        mcols(annot)$isTE <- NULL
    } else
        res <- keep

    if (!returnMask)
        res <- annot[keep]

    res
}

#' @importFrom MatrixGenerics rowRanges
.fetchAnnotationsFromObject <- function(object) {
    annot <- object

    if (is(object, "QuantifyParam"))
        annot <- features(object)
    else if (is(object, "SummarizedExperiment"))
        annot <- rowRanges(object)
    else if (!is(object, "GRangesList") && !is(object, "GRanges"))
        stop(paste("input annotations should be either a 'GenomicRanges',",
                   "'QuantifyParam' or a 'SummarizedExperiment' object."))
    annot
}

#' @importFrom stats setNames
.builDictionary <- function(gr, fuzzy=FALSE) {
  annltr <- gr[gr$repClass=="LTR"]
  # Possible internal id
  intid <- c("int","in", "i")
  intid <- c(paste0("-", intid), paste0("_", intid))
  internal <- grepl(pattern = paste(intid,collapse="|"), annltr$repName, 
                    ignore.case = TRUE)
  ltr <- grepl(pattern = "LTR", annltr$repName, ignore.case = TRUE)
  ltr[internal] <- FALSE # preference goes to internal regions
  # not internal and not ltr: undecided (it could be int or LTR)
  undecid <- !internal & !ltr
  internal[undecid] <- TRUE
  ltr[undecid] <- TRUE
  
  intuniq <- unique(annltr$repName[internal])
  ltruniq <- unique(annltr$repName[ltr])
  
  inside <- character(length(intuniq))
  names(inside) <- intuniq
  outside <- character(length(ltruniq))
  names(outside) <- ltruniq
  
  inout <- .findingEquivalences(ltruniq, intuniq, intid, outside, inside, fuzzy)
  
  namerem <- inout$inside[(names(inout$inside) %in% names(inout$outside))]
  whrem <- !(names(inout$inside) %in% namerem[grep(pattern = "[-_]$",namerem)])
  inout$inside <- inout$inside[whrem]
  whrem <- !(names(inout$outside) %in% namerem[-grep(pattern = "[-_]$",namerem)])
  inout$outside <- inout$outside[whrem]
  
  # When the same LTR is assigned to more than 1 different internal region,
  # the quivalence is added to outside
  whadd <- inout$inside[!(names(inout$inside) %in% inout$outside)]
  inout$outside <- c(inout$outside, setNames(names(whadd), nm = whadd))
  
  # Adding NAs when there is no equivalence for an int region and a LTR
  noteq <- unique(annltr$repName[!(annltr$repName %in% c(inout$inside, inout$outside))])
  noteq <- setNames(rep(NA,length(noteq)), noteq)
  inout$inside <- c(inout$inside, noteq)
  
  inout
}

.findingEquivalences <- function(ltruniq, intuniq, intid, outside, 
                                 inside, fuzzy) {
  # Elements containing "LTR" should have an internal part with a corresponding 
  # I,IN, INT or int suffix
  for (ltr in ltruniq) {
    for (id in c(intid, toupper(intid))) { # c("I","IN", "INT","i","in","int")
      tmpint <- gsub(pattern = "[-_]LTR|LTR", id, ltr)
      if (any(intuniq == tmpint) & ltr != tmpint) {
        inside[tmpint] <- ltr
        outside[ltr] <- tmpint
      }
    }
  }
  # Elements not containing LTR could be matched with a following 
  # "-int", "_int", "_I", "_IN", etc... 
  for (ltr in ltruniq) {
    for (id in c(intid, toupper(intid))) {
      tmpint <- paste0(ltr, id)
      if (any(intuniq == tmpint) & ltr != tmpint) {
        inside[tmpint] <- ltr
        outside[ltr] <- tmpint
      }
    }
  }
  
  ### Combination of the two previous steps
  for (ltr in ltruniq) {
    for (id in c(intid, toupper(intid))) {
      for (idsub in c("I","IN", "INT","i","in","int")) {
        tmpint <- gsub(pattern = "LTR", idsub, ltr)
        tmpint <- paste0(tmpint, id)
        if (any(intuniq == tmpint) & ltr != tmpint) {
          inside[tmpint] <- ltr
          outside[ltr] <- tmpint
        }
      }
    }
  }
  inside <- inside[inside != ""]
  outside <- outside[outside != ""]
  
  # Simplify element names by erasing LTR, I, IN, INT, or symbols like -_, etc...
  elem_liste_ltr <- ltruniq
  ltrpat <- grep(pattern = "^LTR", ltruniq)
  elem_liste_ltr[-ltrpat] <- gsub(pattern = "LTR", ltruniq[-ltrpat], 
                                  replacement = "")
  elem_liste_ltr <- gsub(pattern = "[-_]", replacement = "", elem_liste_ltr)
  names(elem_liste_ltr) <- ltruniq
  elem_liste_int <- gsub(pattern = "int|in", replacement = "", 
                         intuniq, ignore.case = TRUE)
  elem_liste_int <- gsub(pattern = "[-_]I", replacement = "", elem_liste_int, 
                         ignore.case = FALSE)
  elem_liste_int <- gsub(pattern = "[-_]", replacement = "", elem_liste_int)
  names(elem_liste_int) <- intuniq
  whr <- (intuniq %in% ltruniq) & grepl(pattern = "[-_]$", intuniq)
  elem_liste_int <- elem_liste_int[!whr]
  inout <- .simplifyNameSearch(elem_liste_ltr, elem_liste_int, inside, outside)
  
  if (fuzzy) {
    inout <- .getFuzzyEquivalences(inout, elem_liste_int, elem_liste_ltr)
  }
  
  inout
}

.simplifyNameSearch <- function(elem_liste_ltr, elem_liste_int, 
                                inside, outside) {
  ### And we try to match the new names
  el1all <- names(elem_liste_int[is.na(outside[elem_liste_int])])
  el2all <- names(elem_liste_ltr[is.na(outside[elem_liste_ltr]) &
                                   is.na(inside[elem_liste_ltr])])
  for (el1 in el1all) {
    for (el2 in el2all) {
      if (el1 == el2) {
        next
      }
      if (elem_liste_ltr[el2] == elem_liste_int[el1]) {
        if (is.na(inside[el1])) {
          inside[el1] <- el2
        } else {
          inside[el1] <- paste(unique(c(inside[el1], el2)), collapse =":")
        }
        # outside[el2] <- paste(unique(c(inside[el2], el1)), collapse =":")
        if (is.na(outside[el2])) {
          outside[el2] <- el1
        } else {
          outside[el2] <- paste(unique(c(outside[el2], el1)), collapse =":")
        }
      }
    }
  }
  list(inside = inside, outside = outside)
}


.getFuzzyEquivalences <- function(inout, elem_liste_int, elem_liste_ltr) {
  outside <- inout$outside
  inside <- inout$inside
  patterns <- c("[a-z]$","[0-9]$","[a-z0-9][a-z0-9]$")
  for (pat in patterns) {
    for (el1 in sort(names(elem_liste_int))) {
      if (!is.na(outside[el1])) {
        next
      }
      for (el2 in sort(names(elem_liste_ltr))) {
        if ((el1 == el2) | !is.na(outside[el2]) | !is.na(inside[el2])) {
          next
        }
        el1_int <- gsub(pat, "", elem_liste_int[el1], ignore.case = TRUE)
        el2_ltr <- gsub(pat, "", elem_liste_ltr[el2], ignore.case = TRUE)
        if (elem_liste_ltr[el2] == el1_int | elem_liste_int[el1] == el2_ltr) {
          if (is.na(inside[el1])) {
            inside[el1] <- el2
          } else {
            inside[el1] <- paste(unique(c(inside[el1], el2)), collapse =":")
          }
          # outside[el2] <- paste(unique(c(inside[el2], el1)), collapse =":")
          if (is.na(outside[el2])) {
            outside[el2] <- el1
          } else {
            outside[el2] <- paste(unique(c(outside[el2], el1)), collapse =":")
          }
        }
      }
    }
  }
  list(inside = inside, outside = outside)
}


.filterNonTEs <- function(gr) {
  rc <- gr$repClass == "RC" & gr$repFamily == "Helitron"
  gr$repClass[rc] <- "DNA"
  gr$repFamily[rc] <- "RC"
  gr <- gr[gr$repClass %in% c("DNA","SINE","LINE","LTR")]
  gr
}

#' @importFrom GenomicRanges strand
.getConsLength <- function(gr) {
  # Computing length of consensus sequences for each element
  remain <- integer(length(gr))
  remain[as(strand(gr) == "+", "logical")] <- gr$repLeft[as(strand(gr) == "+", "logical")]
  remain[as(strand(gr) == "-", "logical")] <- gr$repStart[as(strand(gr) == "-", "logical")]
  lt <- gr$repEnd + remain
  lt <- split(lt, gr$repName)
  cons_length <- lapply(lt, function(x) as.integer(names(which.max(table(x)))))
  unlist(cons_length)
}

#' @importFrom GenomicRanges seqnames strand
.OCTFTA_parser <- function(chr, gr, cons_length, outside, inside, insert) {
  
  annchrp <- gr[seqnames(gr) == chr & strand(gr) == "+"]
  annchrn <- gr[seqnames(gr) == chr & strand(gr) == "-"]
  if (any(strand(gr) == "*")) { # undefined strand
    annchru <- gr[seqnames(gr) == chr & strand(gr) == "*"]
  }
  
  # Creating GRangesList for the specific chromosome and strand, splitting
  # features based on repName
  annchrp <- split(annchrp, annchrp$repName)
  annchrn <- split(annchrn, annchrn$repName)
  if (any(strand(gr) == "*")) { # undefined strand
    annchru <- split(annchru, annchru$repName)
  }
  
  # checking family is the same for all features with same repName
  # stopifnot(all(unlist(lapply(annchrp, function(gr) length(unique(gr$family)) == 1))))
  # stopifnot(all(unlist(lapply(annchrn, function(gr) length(unique(gr$family)) == 1))))
  
  # Calling .reconstructTEs() to reconstruct TEs
  if (length(annchrp)>0) {
    annchrp_rec <- .reconstructTEs(annchr = annchrp, outside, inside, 
                                   cons_length, insert, minusStrand=FALSE)
  } else {
    annchrp_rec <- annchrp
  }
  if (length(annchrn)>0) {
    annchrn_rec <- .reconstructTEs(annchr = annchrn, outside, inside,
                                   cons_length, insert, minusStrand=TRUE)
  } else {
    annchrn_rec <- annchrn
  }
  annrec <- c(annchrp_rec, annchrn_rec)
  
  if (any(strand(gr) == "*")) {
    if (length(annchru)>0) {
      annchru_rec <- .reconstructTEs(annchr = annchru, outside, inside, 
                                     cons_length, insert, minusStrand=FALSE)
    } else {
      annchru_rec <- annchru
    }
    annrec <- c(annrec, annchru_rec)
  }
  
  annrec
}

#' @importFrom GenomicRanges start mcols "mcols<-"
.reconstructTEs <- function(annchr, outside, inside, cons_length, insert, 
                            minusStrand=FALSE) {
  # For LTRs and internal regions, we check for the presence (in the same 
  # strand) of internal regions and LTRs, respectively, between elements with 
  # the same repName. To do so, we divide features with the same repName
  # (e.g. LTR23) using the start of the corresponding internal region
  # (LTR23-int) as flanking regions.
  
  # Selecting only cases where both the ltr and internal region are in the
  # chromosome
  outsidechr <- outside[(outside %in% names(annchr)) &
                          (names(outside) %in% names(annchr))]
  whint <- which(names(annchr) %in% outsidechr)
  whltr <- which(names(annchr) %in% names(outsidechr))
  # We extract coordinates of int and corresponding ltr (same order)
  annchrint <- annchr[whint]
  annchrltr <- annchr[whltr]
  
  if (length(annchrint) > 0 & length(annchrltr) > 0) {
    splitf2 <- .splitNonContinous(annchrint, annchrltr, outsidechr)
    annchrltrsp <- split(unlist(annchrltr,use.names=FALSE), splitf2$splitfltr2)
    annchrintsp <- split(unlist(annchrint,use.names=FALSE), splitf2$splitfint2)
    annchr <- annchr[-c(whint, whltr)]
    
    # We keep ltr and int separated from the rest of TEs to later perform the
    # reconstruction of full-length or partial ERVs
    
    # Now, features with the same repName are merged together if: they are on
    # the same strand, if the end of the 2nd is closer than 2xconsensus length
    # to the start of the 1st, the start/end of the 2nd is after start of 1st
    # (to avoid merging elements completely contained inside another element),
    # the start/end in the consensus sequence of the 2nd is after the one of 
    # the 1st.
    annchrltrsp2 <- .mergeCloseFeatures(annchrltrsp, cons_length, insert,
                                        minusStrand=minusStrand)
    annchrintsp2 <- .mergeCloseFeatures(annchrintsp, cons_length, insert,
                                        minusStrand=minusStrand)
    
    mcols(annchrltrsp2)$Status <- "LTR"
    mcols(annchrintsp2)$Status <- "int"
  }
  
  # if 'annchrltrsp2' or 'annchrintsp2' are empty, 'annchr2' contains all
  # TEs including LTR and int regions. If 'annchrltrsp2' or 'annchrintsp2' are
  # not empty, still 'annchr' contains LTR and int regions for which no
  # equivalences have been identified (not present in 'outside')
  annchr2 <- .mergeCloseFeatures(annchr, cons_length, insert,
                                 minusStrand=minusStrand)
  
  # Note: All features not identified as LTR or int (from the previously built
  # dictionary with equivalences: 'outside'), are identinfied as "noLTR",
  # however LTRs and int are expected to be inside 'annchr2' since the
  # dictionary 'outside' does not identify all LTR and int regions.
  if (length(annchr2) > 0)
    mcols(annchr2)$Status <- "noLTR"
  
  # ---- Reconstructing full-length and partial ERVs ----
  if (length(annchrint) > 0 & length(annchrltr) > 0) {
    anngrlERVs <- .reconstructERVs(annchrltrsp2, annchrintsp2, outside, inside,
                                   cons_length)
    anngrl <- c(annchr2, anngrlERVs)
  } else {
    anngrl <- annchr2
  }
  opos <- order(min(start(anngrl)), decreasing = FALSE)
  anngrl[opos]
}

#' @importFrom GenomicRanges start mcols "mcols<-"
#' @importFrom Matrix colSums rowSums
.splitNonContinous <- function(annchrint, annchrltr, outsidechr) {
  annchrint2 <- annchrint[outsidechr[names(annchrltr)]]
  annchrltr2 <- split(unlist(annchrltr), 
                      rep(outsidechr[names(annchrltr)], lengths(annchrltr)))
  
  posc_int <- mapply(function(ltrname,intname) 
    outer(start(annchrltr2[[ltrname]]), start(annchrint[[intname]]), "<"),
    ltrname = names(annchrltr2), intname = names(annchrint),
    SIMPLIFY = FALSE)
  
  posc_ltr <- mapply(function(ltrname,intname) 
    outer(start(annchrltr[[ltrname]]), start(annchrint2[[intname]]), "<"),
    ltrname = names(annchrltr), intname = names(annchrint2),
    SIMPLIFY = FALSE)
  
  # 'posc' is a list of matrices in which ltr are rows and int regions are 
  # columns the matrices are TRUE/FALSE depending if the start position of
  # a ltr (row) is smaller than the start position of a int (col)
  # Cases where there is only 1 ltr or internal region in the chr are well
  # addressed
  
  # Creating a factor to split ltr and int according to their position in
  # relation to their equivalent int and ltr, respectively.
  splitfltr <- lapply(posc_ltr, function(x) rowSums(x))
  splitfint <- lapply(posc_int, function(x) colSums(x))
  splitfltr2 <- as.factor(paste(rep(names(splitfltr), lengths(splitfltr)),
                                unlist(splitfltr), sep ="."))
  splitfint2 <- as.factor(paste(rep(names(splitfint), lengths(splitfint)),
                                unlist(splitfint), sep ="."))
  list(splitfltr2 = splitfltr2, splitfint2 = splitfint2)
}

#' @importFrom GenomicRanges GRangesList start width strand end reduce
.mergeCloseFeatures <- function(annchr, cons_length, insert, minusStrand=FALSE) {
  if (insert < 0) {
    # Which elements are separated by less than twice the consensus length?
    annchrwhfirst <- as(lapply(annchr, function(x) x[-1]),"GRangesList")
    stopifnot(identical(names(annchrwhfirst), names(annchr)))
    conslen <- cons_length[do.call("rbind", strsplit(names(annchrwhfirst), 
                                                split = ".", fixed=TRUE))[,1]]
    # yesclose <- (diff(end(annchr)) + width(annchrwhfirst)) < 2*conslen
    yesclose <- (diff(start(annchr)) + width(annchrwhfirst)) < 2*conslen
    
    # Building granges of consensus sequence position
    startcons <- ifelse(strand(unlist(annchr)) == "+",
                        unlist(annchr)$repStart, 
                        unlist(annchr)$repLeft)
    # There are very few cases (for features in "-" strand) where the start is
    # after the end in the consensus sequence. Addressing them by assigning
    # the end position +1 as start 
    disc_start <- startcons > unlist(annchr)$repEnd
    startcons[disc_start] <- unlist(annchr)$repEnd[disc_start] +1
    grcons <- relist(IRanges(start=startcons, 
                             end=unlist(annchr)$repEnd), annchr)
    # Splitting when start/end in chr and consensus sequence of the 1st is not
    # before the 2nd, and elements that are not close enough
    if (minusStrand) {
      yescons <- diff(start(grcons)) < 0 & diff(end(grcons)) < 0
    } else {
      yescons <- diff(start(grcons)) > 0 & diff(end(grcons)) > 0
    }
    yesclosepos <- diff(start(annchr)) > 0 & diff(end(annchr)) > 0 &
      yescons & yesclose
    aggf <- integer(sum(lengths(annchr, use.names = FALSE)))
    whfirst <- c(1, cumsum(lengths(annchr, use.names=FALSE))+1)
    aggf[c(seq_along(aggf))[-whfirst]] <- unlist(cumsum(!yesclosepos))
    annchr2 <- split(unlist(annchr, use.names = FALSE), 
                     paste(rep(names(annchr), lengths(annchr)),aggf, sep ="."))
    annchr2 <- sort(annchr2)
    
    # Note that this implementation does not take into account the start 
    # position of the 1st merged element to see if the distance between the
    # end of a downstream element and the start of an upstream element
    # is smaller than twice the consensus length. Instead, here we take the
    # start of the previous element to compute this distance.
    
  } else {
    insert <- insert + 1 # to adapt to reduce() behaviour
    sp2 <- reduce(annchr, with.revmap=TRUE, min.gapwidth=insert)
    annchr2 <- relist(unlist(annchr, use.names = FALSE), 
                      mcols(unlist(sp2))$revmap)
    names(annchr2) <- paste(rep(names(sp2), lengths(sp2)), 
                            seq_along(annchr2), sep =".")
  }
  annchr2
}

#' @importFrom GenomicRanges start end mcols "mcols<-"
#' @importFrom S4Vectors pc
.reconstructERVs <- function(annchrltrsp2, annchrintsp2, outside, inside,
                             cons_length) {
  ltrtorec <- do.call("rbind", strsplit(names(annchrltrsp2), 
                                        ".", fixed=TRUE))[,1] %in% names(outside)
  inttorec<- do.call("rbind", strsplit(names(annchrintsp2), 
                                       ".", fixed=TRUE))[,1] %in% outside
  # GRangesList with both int and equivalent LTRs
  annchrltrint <- c(annchrltrsp2[ltrtorec], annchrintsp2[inttorec])
  o <- order(min(start(annchrltrint)), decreasing = FALSE)
  annchrltrint <- annchrltrint[o]
  whint <- which(mcols(annchrltrint)$Status == "int")
  namef <- do.call("rbind", strsplit(names(annchrltrint), ".", fixed =TRUE))[,1]
  int <- namef[whint]
  # ltr <- inside[int]
  ltr <- strsplit(inside[int], ":")
  whintup <- whint-1
  whintdo <- whint+1
  # Addressing for when first or last feature in 'annchrltrint' is not an LTR
  # by setting as upstream/downatream feature the int feature itself, causing
  # yesltrup/yesltrdown to be FALSE
  whintup[whintup < 1] <- 1
  whintdo[whintdo > length(annchrltrint)] <- length(annchrltrint)
  # yesltrup <- namef[whintup] == ltr
  # yesltrdown <- namef[whintdo] == ltr
  yesltrup <- .checkEquivalentLTR(namef[whintup], ltr)
  yesltrdown <- .checkEquivalentLTR(namef[whintdo], ltr)
  yesdistup <- (min(start(annchrltrint[whint])) -
                max(end(annchrltrint[whintup]))) < cons_length[namef[whintup]]/2
  # original code uses consensus length of 1st LTR
  yesdistdown <- (min(start(annchrltrint[whintdo])) -
                  max(end(annchrltrint[whint]))) < cons_length[namef[whintup]]/2
  
  # Addressing cases where an LTR is simultaneously a downstream LTR of int1
  # and an upstream LTR of int2. Preference is given to int1.
  whup <- yesltrup & yesdistup
  whdown <- yesltrdown & yesdistdown
  if (sum(whup) > 0 & sum(whdown) > 0) {
    conflictltr <- names(annchrltrint[whintup][whup]) %in%
      names(annchrltrint[whintdo][whdown])
    yesltrup[whup][conflictltr] <- FALSE
    yesdistup[whup][conflictltr] <- FALSE
    whup <- yesltrup & yesdistup
    whdown <- yesltrdown & yesdistdown
  }
  fulllength <- whup & whdown
  partial <- whup | whdown
  partial[fulllength] <- FALSE
  
  anngrlERVs <- .getFulllength_Partial_ERVs(fulllength, partial, annchrltrint, 
                                            whintup, whint, whintdo, yesltrup,
                                            yesdistup, yesltrdown, yesdistdown,
                                            annchrltrsp2, annchrintsp2, ltrtorec,
                                            inttorec)
  anngrlERVs
}

#' @importFrom GenomicRanges start end mcols "mcols<-"
#' @importFrom S4Vectors pc
.getFulllength_Partial_ERVs <- function(fl, pt, annchrltrint, whintup,
                                        whint, whintdo, yesltrup, yesdistup,
                                        yesltrdown, yesdistdown, annchrltrsp2,
                                        annchrintsp2, ltrtorec, inttorec) {
  # Creating GRangesList with full-length ERVs
  if (any(fl)) {
    fulllength_grl <- pc(annchrltrint[whintup][fl], annchrltrint[whint][fl],
                         annchrltrint[whintdo][fl])
    mcols(fulllength_grl)$Status <- "full-lengthLTR"
    # assigning names of int element
    names(fulllength_grl) <- names(annchrltrint[whint][fl])
    # Making sure no partially reconstructed ERV contains fragments from a 
    # full-length ERV
    disc <- names(annchrltrint[whintup][pt]) %in% names(annchrltrint[whintdo][fl])
    disc2 <- names(annchrltrint[whintdo][pt]) %in% names(annchrltrint[whintup][fl])
    pt[pt][disc | disc2] <- FALSE
  } else {
    fulllength_grl <- GRangesList()
  }
  ptup <- pt & yesltrup & yesdistup
  ptdown <- pt & yesltrdown & yesdistdown
  
  # Creating GRangesList of partial ERVs
  if (any(ptup)) {
    partial_grl1 <- pc(annchrltrint[whintup][ptup], annchrltrint[whint][ptup])
    mcols(partial_grl1)$Status <- "partialLTR_up"
    # assigning names of int element
    names(partial_grl1) <- names(annchrltrint[whint][ptup])
    # removing LTR for reconstruction with an upstream int if it has been 
    # previously reconstructed with a downstream int
    disc <- names(annchrltrint[whintdo][ptdown]) %in%
      names(annchrltrint[whintup][ptup])
    pt[ptdown][disc] <- FALSE
    ptdown <- pt & yesltrdown & yesdistdown
  } else {
    partial_grl1 <- GRangesList()
  }
  if (any(ptdown)) {
    partial_grl2 <- pc(annchrltrint[whint][ptdown],annchrltrint[whintdo][ptdown])
    mcols(partial_grl2)$Status <- "partialLTR_down"
    # assigning names of int element
    names(partial_grl2) <- names(annchrltrint[whint][ptdown])
  } else {
    partial_grl2 <- GRangesList()
  }
  # names() of partial_grl2 are already int elements
  partial_grl <- c(partial_grl1, partial_grl2)
  
  # Which are not solo int or LTR (have been reconstructed)?
  whyesrec <- c(whint[fl], c(whintdo)[fl], c(whintup)[fl], whint[ptup], 
                whint[ptdown], c(whintup)[ptup], c(whintdo)[ptdown])
  
  # GRangesList with reconstructed ERVs and solo LTR and int
  anngrlERVs <- c(annchrltrsp2[!ltrtorec], annchrintsp2[!inttorec],
                  fulllength_grl, partial_grl, annchrltrint[-whyesrec])
  anngrlERVs
}

.checkEquivalentLTR <- function(namef, ltr) {
  yesequiv <- logical()
  for (i in seq_along(namef)) {
    yesequiv <- c(yesequiv, namef[i] %in% ltr[[i]])
  }
  yesequiv
}

#' @importFrom GenomicRanges width mcols "mcols<-"
.getRelLength <- function(nTEs, cons_length, inout, annrec) {
  whERVs <- mcols(annrec)$Status %in% c("partialLTR_up", "partialLTR_down",
                                        "full-lengthLTR", "int")
  ltrconsl <- cons_length[names(inout$outside)[match(nTEs[whERVs], inout$outside)]]
  cons_length_full <- cons_length[nTEs[whERVs]] + 2*ltrconsl
  # cons_length_full <- cons_length[nTEs[whERVs]] + 2*cons_length[inside[nTEs[whERVs]]]
  rlenERVs <- sum(width(annrec[whERVs])) / cons_length_full
  
  # for solo LTRs and the rest of TEs the consensus length is obtained from the
  # positions in the consensus sequence
  rlenOther <- sum(width(annrec[!whERVs])) / cons_length[nTEs[!whERVs]]
  
  rlen <- numeric(length(annrec))
  rlen[whERVs] <- rlenERVs
  rlen[!whERVs] <- rlenOther
  # addressing relative lengths > 1 (due to assembled fragments that partially
  # overlap between them or consensus lengths shorter than the single fragment)
  rlen[rlen > 1] <- 1
  
  rlen
}


#' atena annotation parser of RepeatMasker annotations
#' @param gr A \link[GenomicRanges:GRanges-class]{GRanges} object with
#' RepeatMasker annotations from \link[AnnotationHub]{AnnotationHub}
#' 
#' @param strict (Default FALSE) A logical; if TRUE, the 80-80 rule is applied,
#' i.e. only copies with more than 80% identity to the reference
#' and more than 80 bp long are reported.
#' 
#' @param insert (Default 1000L) An integer > 0. Fragments are assembled 
#' together if the distance between their closest extremities
#' is equal or less than \code{insert}. When \code{insert} = 0, two fragments
#' are assembled if they are in contact next to each other.
#'
#' @return A \link[GenomicRanges:GRangesList-class]{GRangesList} object.
#'         
#' @details 
#' atena annotation parser of RepeatMasker annotations.
#' Parses RepeatMasker annotations from UCSC by assembling together fragments
#' from the same transposable element (TE) that are close enough (determined
#' by the \code{insert} parameter). For TEs from the LTR class, the parser
#' tries to reconstruct full-length, when possible, or partial TEs following
#' the LTR - internal region - LTR structure. Equivalences between LTR and
#' internal regions are found by, first, identifying LTR regions (those with
#' the "LTR" substring in their name) and internal regions (those with a
#' suffix such as "-int", "-I", etc.). Then, LTR are assigned to internal
#' regions for which the comparison of the two names are has a higher number
#' of equal consecutive characters.
#'
#' @examples
#' rmskat <- annotaTEs(genome="dm6", parsefun=rmskatenaparser,
#'                     strict=FALSE)
#' 
#' @aliases rmskatenaparser
#' @rdname rmskatenaparser
#' @name rmskatenaparser
#' @importFrom GenomicRanges strand width mcols "mcols<-" seqnames
#' @importFrom IRanges mean
#' @importFrom GenomeInfoDb seqlevels
#' @export
rmskatenaparser <- function(gr, strict= FALSE, insert=1000) {
  if (!is(gr, "GRanges"))
    stop("'gr' should be a GRanges object with RepeatMasker annotations")
  
  if (length(gr) == 0)
    stop("'gr' is empty")
  
  if (!is.integer(insert) & !is.numeric(insert))
    stop("'insert' must be an integer value")
  
  inout <- .builDictionary_at(gr)
  
  gr <- .filterNonTEs(gr)
  gr$repStart <- abs(gr$repStart)
  gr$repLeft <- abs(gr$repLeft)
  cons_length <- .getConsLength(gr)
  
  fsplit <- paste(mcols(gr)$repName, seqnames(gr), strand(gr), sep = ";")
  ann <- split(gr, fsplit)
  mcols(ann)$repName <- do.call("rbind", strsplit(names(ann), split = ";", 
                                                  fixed =TRUE))[,1]
  annrec <- .reconstructTEs_at(ann = ann, inout$outside, inout$inside,
                               cons_length, insert)
  # -- Filtering based on 'strict' option
  if (strict) {
    maxdiv <- mean(relist(mcols(unlist(annrec))$milliDiv, annrec)) < 200
    tokeep <- sum(width(annrec)) >= 80 & maxdiv
    annrec <- annrec[tokeep]
  }
  
  # -- Computation relative length --
  # for full-length and partial ERVs, and ERVs with only internal region, the
  # consensus length is considered that of the full-length TE: LTR + int + LTR
  nTEs <- do.call("rbind", strsplit(names(annrec), ";", fixed=TRUE))[,1]
  mcols(annrec)$RelLength <- .getRelLength(nTEs, cons_length, inout, annrec)
  
  # -- Simplifying TE name: subfamily name + number of TE separated by "." --
  names(annrec) <- paste(nTEs, seq_along(annrec), sep = ".")
  
  # -- Including TE class name in GRangesList 
  mcols(annrec)$Class <- unlist(unique(relist(unlist(annrec)$repClass, annrec)))
  
  annrec
}

#' @importFrom stats aggregate
#' @importFrom Matrix colSums
.builDictionary_at <- function(gr, thmatch=5, threl = 0.5) {
  
  # identifying LTR and int regions
  annltr <- gr[gr$repClass=="LTR"]
  # Possible internal id
  intid <- c("int","in", "i")
  intid <- c(paste0("-", intid), paste0("_", intid))
  internal <- grepl(pattern = paste(intid,collapse="|"), annltr$repName, 
                    ignore.case = TRUE)
  ltr <- grepl(pattern = "LTR", annltr$repName, ignore.case = TRUE)
  ltr[internal] <- FALSE # preference goes to internal regions
  # not internal and not ltr: undecided (it could be int or LTR)
  # here we discard undecided elements
  undecid <- !internal & !ltr
  internal[undecid] <- FALSE
  ltr[undecid] <- FALSE
  
  intuniq <- unique(annltr$repName[internal])
  ltruniq <- unique(annltr$repName[ltr])
  
  inside <- character(length(intuniq))
  names(inside) <- intuniq
  outside <- character(length(ltruniq))
  names(outside) <- ltruniq
  
  mm <- NULL
  for (ltrn in names(outside)) {
    nmatch <- vapply(names(inside), function(intn) {
      max(nchar(abbreviate(c(intn, ltrn), minlength = 1, use.classes=FALSE)))
    }, FUN.VALUE = integer(1L))
    nmatch[nmatch < c(thmatch +1)] <- 0
    nmatch <- nmatch / nchar(ltrn)
    mm <- cbind(mm, nmatch)
  }
  colnames(mm) <- names(outside)
  rownames(mm) <- names(inside)
  # at least 50% of characters in LTR name have to consecutively match the
  # internal region name
  mm[mm < threl] <- 0L
  whmin <- colSums(mm) > 0
  outside[whmin] <- names(inside)[apply(mm[,whmin],2, which.max)]
  aginside <- aggregate(names(outside[outside != ""]), 
                        by = list(outside[outside != ""]),
                        FUN= paste, collapse = ":")
  inside[aginside$Group.1] <- aginside$x
  outside <- outside[outside != ""]
  
  list(inside = inside, outside = outside)
}

#' @importFrom GenomicRanges start mcols "mcols<-"
#' @importFrom GenomeInfoDb seqlevels
.reconstructTEs_at <- function(ann, outside, inside, cons_length, insert) {
  # For LTRs and internal regions, we check for the presence (in the same 
  # strand) of internal regions and LTRs, respectively, between elements with 
  # the same repName. To do so, we divide features with the same repName
  # (e.g. LTR23) using the start of the corresponding internal region
  # (LTR23-int) as flanking regions.
  
  # Selecting only cases where both the ltr and internal region are in the
  # chromosome
  stchr <- c(paste(seqlevels(ann), c("+"), sep =";"),
             paste(seqlevels(ann), c("-"), sep =";"))
  outside_stchr <- paste(rep(outside, each = length(stchr)), 
                         rep(stchr, length(outside)), sep =";")
  names(outside_stchr) <-  paste(rep(names(outside), each = length(stchr)), 
                                 rep(stchr, length(outside)), sep =";")
  outsidechr <- outside_stchr[(outside_stchr %in% names(ann)) &
                                (names(outside_stchr) %in% names(ann))]
  whint <- which(names(ann) %in% outsidechr)
  whltr <- which(names(ann) %in% names(outsidechr))
  # We extract coordinates of int and corresponding ltr (same order)
  annint <- ann[whint]
  annltr <- ann[whltr]
  
  if (length(annint) > 0 & length(annltr) > 0) {
    splitf2 <- .splitNonContinous_at(annint, annltr, outsidechr)
    annltrsp <- split(unlist(annltr, use.names = FALSE), splitf2$splitfltr2)
    annintsp <- split(unlist(annint, use.names = FALSE), splitf2$splitfint2)
    ann <- ann[-c(whint, whltr)]
    # We keep ltr and int separated from the rest of TEs to later perform the
    # reconstruction of full-length or partial ERVs
    
    # Now, features with the same repName, strand and chromosome are merged
    # together if they are close enough, according to the 'insert' parameter
    annltrsp2 <- .mergeCloseFeatures_at(annltrsp, cons_length, insert)
    annintsp2 <- .mergeCloseFeatures_at(annintsp, cons_length, insert)
    mcols(annltrsp2)$Status <- "LTR"
    mcols(annintsp2)$Status <- "int"
  }
  ann2 <- .mergeCloseFeatures_at(ann, cons_length, insert)
  
  # Note: All features not identified as LTR or int (from the previously built
  # dictionary with equivalences: 'outside'), are identinfied as "noLTR",
  # however LTRs and int are expected to be inside 'annchr2' since the
  # dictionary 'outside' does not identify all LTR and int regions.
  mcols(ann2)$Status <- "noLTR"
  
  # ---- Reconstructing full-length and partial ERVs ----
  if (length(annint) > 0 & length(annltr) > 0) {
    anngrlERVs <- .reconstructERVs_at(annltrsp2, annintsp2, outsidechr,
                                      outside, inside, cons_length)
    anngrl <- c(ann2, anngrlERVs)
  } else {
    anngrl <- ann2
  }
  chr <- do.call("rbind", strsplit(names(anngrl), ";", fixed=TRUE))[, 2]
  opos <- order(chr, min(start(anngrl)), decreasing=FALSE)
  anngrl[opos]
}

#' @importFrom GenomicRanges start mcols "mcols<-"
#' @importFrom Matrix colSums rowSums
.splitNonContinous_at <- function(annint, annltr, outsidechr) {
  annint2 <- annint[outsidechr[names(annltr)]]
  annltr2 <- split(unlist(annltr), 
                   rep(outsidechr[names(annltr)], lengths(annltr)))
  
  posc_int <- mapply(function(ltrname, intname) 
    outer(start(annltr2[[ltrname]]), start(annint[[intname]]), "<"),
    ltrname = names(annltr2), intname = names(annint),
    SIMPLIFY = FALSE)
  
  posc_ltr <- mapply(function(ltrname, intname) 
    outer(start(annltr[[ltrname]]), start(annint2[[intname]]), "<"),
    ltrname = names(annltr), intname = names(annint2),
    SIMPLIFY = FALSE)
  
  # 'posc' is a list of matrices in which ltr are rows and int regions are 
  # columns the cells are TRUE/FALSE depending if the start position of
  # a ltr (row) is smaller than the start position of a int (col)
  
  # Creating a factor to split ltr and int according to their position in
  # relation to their equivalent int and ltr, respectively.
  splitfltr <- lapply(posc_ltr, function(x) rowSums(x))
  splitfint <- lapply(posc_int, function(x) colSums(x))
  splitfltr2 <- as.factor(paste(rep(names(splitfltr), lengths(splitfltr)),
                                unlist(splitfltr), sep = "."))
  splitfint2 <- as.factor(paste(rep(names(splitfint), lengths(splitfint)),
                                unlist(splitfint), sep = "."))
  list(splitfltr2 = splitfltr2, splitfint2 = splitfint2)
}

#' @importFrom GenomicRanges GRangesList start width strand end reduce
.mergeCloseFeatures_at <- function(annsp, cons_length, insert) {
  insert <- insert + 1 # to adapt to reduce() behaviour
  sp2 <- reduce(annsp, with.revmap=TRUE, min.gapwidth=insert)
  annsp2 <- relist(unlist(annsp, use.names=FALSE), mcols(unlist(sp2))$revmap)
  names(annsp2) <- paste(rep(names(sp2), lengths(sp2)), 
                         seq_along(annsp2), sep =".")
  annsp2
}

#' @importFrom GenomicRanges start end mcols "mcols<-"
#' @importFrom S4Vectors pc
.reconstructERVs_at <- function(annltrsp2, annintsp2, outsidechr, 
                                outside, inside, cons_length) {
  ltrtorec <- do.call("rbind", strsplit(names(annltrsp2), 
                                        ".", fixed=TRUE))[,1] %in% names(outsidechr)
  inttorec<- do.call("rbind", strsplit(names(annintsp2), 
                                       ".", fixed=TRUE))[,1] %in% outsidechr
  # GRangesList with both int and equivalent LTRs
  annltrint <- c(annltrsp2[ltrtorec], annintsp2[inttorec])
  
  # Setting order based on strand, chr and start position of element 
  chr <- do.call("rbind", strsplit(names(annltrint), ";", fixed=TRUE))[,2]
  st <- do.call("rbind", strsplit(names(annltrint), ";", fixed=TRUE))[,3]
  st <- do.call("rbind", strsplit(st, ".", fixed=TRUE))[,1]
  o <- order(chr, st, min(start(annltrint)), decreasing = FALSE)
  annltrint <- annltrint[o]
  whint <- which(mcols(annltrint)$Status == "int")
  namef <- do.call("rbind", strsplit(names(annltrint), ".", fixed =TRUE))[,1]
  int <- namef[whint]
  ltr <- split(names(outsidechr), f = outsidechr)
  ltr <- ltr[int]
  whintup <- whint-1
  whintdo <- whint+1
  # Addressing for when first or last feature in 'annltrint' is not an LTR
  # by setting as upstream/downatream feature the int feature itself, causing
  # yesltrup/yesltrdown to be FALSE
  whintup[whintup < 1] <- 1
  whintdo[whintdo > length(annltrint)] <- length(annltrint)
  yesltrup <- .checkEquivalentLTR(namef[whintup], ltr)
  yesltrdown <- .checkEquivalentLTR(namef[whintdo], ltr)
  namef_simp <- do.call("rbind", strsplit(namef, ";", fixed =TRUE))[,1]
  yesdistup <- (min(start(annltrint[whint])) - max(end(annltrint[whintup]))) <
    cons_length[namef_simp[whintup]]/2
  yesdistdown <- (min(start(annltrint[whintdo])) - max(end(annltrint[whint]))) <
    cons_length[namef_simp[whintup]]/2
  
  # Addressing cases where an LTR is simultaneously a downstream LTR of int1
  # and an upstream LTR of int2. Preference is given to int1.
  whup <- yesltrup & yesdistup
  whdown <- yesltrdown & yesdistdown
  if (sum(whup) > 0 & sum(whdown) > 0) {
    conflictltr <- names(annltrint[whintup][whup]) %in%
      names(annltrint[whintdo][whdown])
    yesltrup[whup][conflictltr] <- FALSE
    yesdistup[whup][conflictltr] <- FALSE
    whup <- yesltrup & yesdistup
    whdown <- yesltrdown & yesdistdown
  }
  fulllength <- whup & whdown
  partial <- whup | whdown
  partial[fulllength] <- FALSE
  
  anngrlERVs <- .getFulllength_Partial_ERVs(fulllength, partial, annltrint,
                                            whintup, whint, whintdo, yesltrup,
                                            yesdistup, yesltrdown, yesdistdown,
                                            annltrsp2, annintsp2, ltrtorec,
                                            inttorec)
  anngrlERVs
}

