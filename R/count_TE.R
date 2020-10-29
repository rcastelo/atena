#' Quantification of Transposable Elements based on ERVmap pipeline
#'
#' \code{count_TE} quantifies the expression of transposable elements (TE) at
#' the locus level. The approach used by this function is based on the pipeline
#' described by
#' \href{https://doi.org/10.1073/pnas.1814589115}{Tokuyama et al. (2018)},
#' in which reads are processed following
#' conservative filtering criteria to provide reliable raw count data for each
#' genomic region encompassing a TE. Only those alignments that pass the filters
#' are considered for TE expression quantification.
#'
#' Briefly, alignments need to pass three filtering criteria in order to be
#' included in the quantification of TEs. These three filters are:
#' \enumerate{
#'   \item Sum of hard and soft clipping bases to the length of the read must
#'   be lower than 0.02.
#'   \item The ratio of the edit distance to the sequence read length must be
#'   less than 0.02.
#'   \item The difference between the alignment score (AS tag) and the
#'   suboptimal alignment score (XS tag) must be equal or higher than 5.
#' }
#'
#' As mentioned, this function requires the tag containing the suboptimal
#' alignment score in order to correctly perform the filtering. This is why,
#' for those datasets in which this information is not provided in the binary
#' 'BAM' file, the \code{\link{getSuboptimalAlignScore}} function is called
#' internally to compute the suboptimal alignment scores.
#'
#' Once the filtering criteria is applied to the alignments, only those
#' meeting the criteria are considered to obtain the read counts for each
#' element in the annotation file. The
#' \code{\link[GenomicAlignments]{summarizeOverlaps}} function is used for this
#' purpose.
#'
#' @param path_bam_files_sorted Character vector indicating the location where
#'   the sorted by Query-name BAM files resulting from the \code{\link{sort_bam}}
#'   function are located.
#'
#' @param strandMode (Default 1) Numeric vector which can take values 0, 1 or 2.
#'   The strand mode is a per-object switch on
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   objects that controls the behavior of the strand getter. See
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   class for further detail. If \code{singleEnd = TRUE}, then use either
#'   \code{strandMode = NULL} or do not specify the \code{strandMode} parameter.
#'
#' @param ann_file (Default NULL) A single string containing the path and name
#'   of the endogenous retrovirus (ERVs) or TE annotation file (including the
#'   extension of the file). The allowed file formats are BED and GTF / GFF
#'   version 2. When \code{ann_file} is not specified (\code{ann_file = NULL}),
#'   the human endogenous retroviruses annotations published by
#'   \href{https://doi.org/10.1073/pnas.1814589115}{Tokuyama et al. (2018)}
#'   will be used instead. These annotations are not strand-specific, so
#'   \code{ignore_strand} must be \code{ignore_strand = TRUE}.
#'
#' @param eval_unique_r (Default TRUE) Logical value indicating whether to apply
#' the alignment filters to unique reads (TRUE) or not (FALSE). These filters,
#' which are always applied to multi-mapping reads, are optional for unique
#' reads. If TRUE, the unique reads not passing one or more filters from the
#' ERVmap pipeline, except for the "AS - XS >= 5" filter, will be discarded to
#' compute TEs expression.
#'
#' @param ignore_strand (Default TRUE) A logical which defines if the strand
#' should be taken into consideration when computing the overlap between reads
#' and TEs/ERVs in the annotations. When \code{ignore_strand = FALSE}, the
#' \code{\link[GenomicAlignments]{summarizeOverlaps}} function will only
#' consider those reads selected after filtering which overlap the TE or
#' ERV on the same strand. On the contrary, when \code{ignore_strand = TRUE},
#' the \code{\link[GenomicAlignments]{summarizeOverlaps}} function will count
#' any alignment which overlaps with the element in the annotations regardless
#' of the strand. For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}}.
#'
#' @param singleEnd (Default FALSE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @return A `SummarizedExperiment` object containing the read counts for
#' each sample and element in the annotations file. Each sample is
#' represented as a column and each row corresponds to an element in the
#' annotations. Column names are extracted from the binary 'BAM' file name of
#' each sample.
#'
#' @seealso
#' \code{\link[GenomicAlignments]{GAlignmentPairs-class}},
#' \code{\link[GenomicAlignments]{summarizeOverlaps}},
#' \code{\link{getSuboptimalAlignScore}}
#'
#' @examples
#' path_bam_files_sorted <- system.file("extdata", package = "atena")
#'
#' jong_HERV_se <- count_TE(path_bam_files_sorted, eval_unique_r = TRUE,
#'                          ignore_strand = TRUE), singleEnd = TRUE)
#'
#' @references
#' Tokuyama M et al. ERVmap analysis reveals genome-wide transcription of human
#' endogenous retroviruses. PNAS. 2018;115(50):12565-12572. DOI:
#' \url{https://doi.org/10.1073/pnas.1814589115}
#'
#' @export
#' @import GenomicAlignments
#' @import Rsamtools
#' @import rtracklayer
#' @import S4Vectors
#' @import SummarizedExperiment



count_TE <- function(path_bam_files_sorted, singleEnd = FALSE, strandMode = 1,
                     ann_file = NULL, eval_unique_r = TRUE,
                     ignore_strand = TRUE) {

  # -- Importing HERVs coordinates

  if (is.null(ann_file)) {
    ERV_ann <- ERVmap_ann

  } else {

    # Obtaining the format of the ERVs annotations
    ERV_format <- tail(strsplit(ann_file, split = ".", fixed = TRUE)[[1]],
                       n = 1)

    # Reading the annotations file according to the format
    if (ERV_format == "bed") {
      ERV_ann <- rtracklayer::import.bed(ann_file)

    } else if (ERV_format == "gtf" | ERV_format == "gff") {
      ERV_ann <- rtracklayer::readGFF(ann_file)
      ERV_ann_granges <- GenomicRanges::makeGRangesFromDataFrame(ERV_ann)
      mcols(ERV_ann_granges) <- ERV_ann$gene_id
      colnames(mcols(ERV_ann_granges))[1] <- "name"
      ERV_ann <- ERV_ann_granges
    }

    # Changing "seqnames" to match the ones from the reads
    if (!all(grepl(x = levels(seqnames(ERV_ann)),
                   pattern = "chr", fixed = TRUE))) {
      GenomeInfoDb::seqlevels(ERV_ann) <- levels(as.factor(paste("chr",
                                                                 seqnames(ERV_ann),
                                                   sep = '')))
    }
  }


  # --- Single-end analysis --- #
  if (singleEnd) {
    # -- Reading some alignments of the BAM files to see if they contain the
    # "NH", "XS", "AS" and "NM" tags
    bam_files <- list.files(path_bam_files_sorted, pattern = "\\.bam$",
                            full.names = TRUE)
    bfl <- BamFileList(bam_files, yieldSize = 100)

    sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                           isDuplicate = FALSE,
                           isNotPassingQualityControls = FALSE)

    param <- ScanBamParam(flag = sbflags,
                          tag = c("nM", "NM", "AS", "NH", "XS"))

    NH_tag <- NULL
    NM_tag <- NULL
    nM_tag <- NULL
    AS_tag <- NULL
    XS_tag <- NULL

    for (i in 1:length(bfl)) {

      file <- bfl[[i]]
      open(file)
      r_test <- readGAlignments(file, param = param)
      close(file)

      # Testing, for each sample, if the files contain the different tags
      NH_tag[i] <- ifelse(all(is.na(mcols(r_test)$NH)), FALSE, TRUE)
      NM_tag[i] <- ifelse(all(is.na(mcols(r_test)$NM)), FALSE, TRUE)
      nM_tag[i] <- ifelse(all(is.na(mcols(r_test)$nM)), FALSE, TRUE)
      AS_tag[i] <- ifelse(all(is.na(mcols(r_test)$AS)), FALSE, TRUE)
      XS_tag[i] <- ifelse(all(grepl("\\d", (mcols(r_test)$XS))), TRUE, FALSE)

    }

    # -- Setting the parameters for reading the BAM files
    bfl <- BamFileList(bam_files, yieldSize = 100000)

    sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                           isDuplicate = FALSE,
                           isNotPassingQualityControls = FALSE,
                           isSupplementaryAlignment = FALSE)

    # If neither NH nor XS tags are available, the suboptimal alignment score
    # cannot be computed and the 3rd filter cannot be applied. Plus, unique
    # reads cannot be differentiated from multi-mapping reads
    if (!all(NH_tag) & !all(XS_tag)) {

      if (eval_unique_r == FALSE) {
        stop("Error in 'eval_unique_r = FALSE': Neither the NH tag nor the XS tag (from BWA) are provided in the BAM files. Unique reads cannot be differentiated from multi-mapping reads, therefore unique reads must also be filtered.")
      }

      warning("Neither the NH tag nor the XS tag (from BWA) are provided in the BAM files. At least one of the two tags is needed to obtain the suboptimal alignment score. The 3rd filter (based on the suboptimal alignment score) will not be applied.")
      param <- ScanBamParam(flag = sbflags,
                            tag = c("nM", "NM", "AS", "NH", "XS"))

    } else if (all(NH_tag) & !all(XS_tag)) {
      param <- ScanBamParam(flag = sbflags,
                            tag = c("nM", "NM", "AS", "NH", "XS"),
                            tagFilter = list(NH = c(2:10000)))

    } else if (!all(NH_tag) & all(XS_tag)) {
      param <- ScanBamParam(flag = sbflags,
                            tag = c("nM", "NM", "AS", "NH", "XS"),
                            tagFilter = list(XS = c(1:10000)))

    } else {
      stop("Error: an error occurred when reading the tags from the BAM files.")
    }

    table_counts <- matrix(nrow = length(ERV_ann), ncol = length(bfl))
    colnames(table_counts) <- gsub("\\.bam", "", names(bfl))


    for (i in 1:length(bfl)) {

      file <- bfl[[i]]
      open(file)

      while (length(r <- readGAlignments(file, param = param))) {

        # Compute the suboptimal alignment score (equivalent to XS tag in BWA)
        # only if the score is not available in the BAM file and the NH tag is
        # available
        if (all(NH_tag) & !all(XS_tag)) {

          if (is.null(names(r))) {
            stop("The GAlignmentPairs object must have names")
          }

          r <- atena::getSuboptimalAlignScore(r, strandMode = strandMode)
        }


        # The first filter is always applied
        cigar <- explodeCigarOpLengths(cigar(r), ops = c("H","S"))
        SH_clipping <- lapply(cigar, function(x) sum(unlist(x)))

        # The second filter is always applied but when the NM tag is not
        # available the nM tag can be used instead.
        if (all(NM_tag)) {
          NM_filter <- (mcols(r)$NM / qwidth(r)) < 0.02
        } else if (all(nM_tag)) {
          NM_filter <- (mcols(r)$nM / qwidth(r)) < 0.02
        } else {
          stop("Error: Neither the NM tag nor the nM tag are available in the BAM file. At least one of them is needed to apply the 2nd filter.")
        }

        # The third filter (AS - XS >= 5) cannot be applied when neither the NH
        # tag nor the XS tag are available
        if (!all(NH_tag) & !all(XS_tag)) {
          AS_XS_filter <- rep(TRUE, length(r))
        } else {
          AS_XS_filter <- (mcols(r)$AS - mcols(r)$XS) >= 5
          AS_XS_filter[is.na(AS_XS_filter)] <- FALSE
        }

        # Filtering alignments considering all filters
        r_to_keep <- ((unlist(SH_clipping) / qwidth(r)) < 0.02) &
          NM_filter & AS_XS_filter

        r_total <- r[r_to_keep]
        rm(r)

        overlap <- summarizeOverlaps(features = ERV_ann, reads = r_total,
                                     ignore.strand = ignore_strand,
                                     mode = Union, inter.feature = FALSE,
                                     singleEnd = TRUE)
        rm(r_total)

        if (all(is.na(table_counts[,i]))) {
          table_counts[,i] <- assay(overlap)
        } else {
          table_counts[,i] <- table_counts[,i] + assay(overlap)
        }


      }

      close(file)

      if (any(c(all(NH_tag),all(XS_tag)))) {
        ## If the NH and XS are not provided, in the previous step both unique
        ## reads and multimapping reads have already been filtered
        ## If not, if the user wants to filter unique reads, they will be
        ## filtered (only the first 2 filters). Instead, if the user does not
        ## want to filter unique reads, they will not be filtered, but they
        ## will be considered for computing HERV raw counts

        # Now the alignments corresponding to unique reads are read
        sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                               isSupplementaryAlignment = FALSE,
                               isDuplicate = FALSE,
                               isNotPassingQualityControls = FALSE)

        if (!all(NH_tag) & all(XS_tag)) {

          # If the NH is not available but the XS is, unique reads can be
          # read by reading those in which XS = 0.
          param <- ScanBamParam(flag = sbflags,
                                tag = c("nM", "NM", "AS", "NH", "XS"),
                                tagFilter = list(XS = 0))
        } else {
          param <- ScanBamParam(flag = sbflags,
                                tag = c("nM","NM","AS","NH","XS"),
                                tagFilter = list(NH = 1)) #Read only unique reads
        }

        open(file)

        while (length(r <- readGAlignments(file, param = param))) {

          # If the user wants the unique reads to pass the filters in ERVmap,
          # run ERVmap (except for the "AS - XS > 5" filter) on the unique reads
          if (eval_unique_r) {

            cigar <- explodeCigarOpLengths(cigar(r), ops = c("H","S"))
            SH_clipping <- lapply(cigar, function(x) sum(unlist(x)))

            if (all(NM_tag)) {
              NM_filter <- (mcols(r)$NM / qwidth(r)) < 0.02
            } else if (all(nM_tag)) {
              NM_filter <- (mcols(r)$nM / qwidth(r)) < 0.02
            } else {
              stop("Error: Neither the NM tag nor the nM tag are available. At least one of them is needed to apply the 2nd filter.")
            }

            # Filtering alignments considering all filters
            r_to_keep <- ((unlist(SH_clipping) / qwidth(r)) < 0.02) & NM_filter
            r_unique <- r[r_to_keep]
            rm(r)

            overlap <- summarizeOverlaps(features = ERV_ann, reads = r_unique,
                                         ignore.strand = ignore_strand,
                                         mode = Union, inter.feature = FALSE,
                                         singleEnd = TRUE)
            rm(r_unique)

            # If the user does not want the unique reads to be filtered, but
            # instead wants to include all unique alignments in the analysis
            # (discarding unmapped, unpaired, duplicated and not passing
            # quality control reads), then:
          } else {
            overlap <- summarizeOverlaps(features = ERV_ann, reads = r,
                                         ignore.strand = ignore_strand,
                                         mode = Union, inter.feature = FALSE,
                                         singleEnd = TRUE)
            rm(r)
          }

          if (all(is.na(table_counts[,i]))) {
            table_counts[,i] <- assay(overlap)
          } else {
            table_counts[,i] <- table_counts[,i] + assay(overlap)
          }
        }
        gc()
        close(file)

      }

      sample_name <- gsub("\\.bam", "", names(bfl)[i])
      colnames(table_counts)[i] <- sample_name

      message(paste("Sample", sample_name, "done!", sep = " "))
    }


    # --- Paired-end analysis --- #
  } else {
    # -- Reading some alignments of the BAM files to see if they contain the
    # "NH", "XS", "AS" and "NM" tags
    bam_files <- list.files(path_bam_files_sorted, pattern = "\\.bam$",
                            full.names = TRUE)
    bfl <- BamFileList(bam_files, asMates = TRUE, yieldSize = 100)

    sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                           isProperPair = TRUE,
                           isDuplicate = FALSE,
                           isSupplementaryAlignment = FALSE,
                           isNotPassingQualityControls = FALSE)

    param <- ScanBamParam(flag = sbflags,
                          tag = c("nM", "NM", "AS", "NH", "XS"))

    NH_tag <- NULL
    NM_tag <- NULL
    nM_tag <- NULL
    AS_tag <- NULL
    XS_tag <- NULL

    for (i in 1:length(bfl)) {

      file <- bfl[[i]]
      open(file)
      r_test <- readGAlignmentPairs(file, use.names = TRUE,
                                    param = param,
                                    strandMode = strandMode)
      close(file)

      # Testing, for each sample, if the files contain the different tags
      NH_tag[i] <- ifelse(all(is.na(mcols(first(r_test))$NH)), FALSE, TRUE)
      NM_tag[i] <- ifelse(all(is.na(mcols(first(r_test))$NM)), FALSE, TRUE)
      nM_tag[i] <- ifelse(all(is.na(mcols(first(r_test))$nM)), FALSE, TRUE)
      AS_tag[i] <- ifelse(all(is.na(mcols(first(r_test))$AS)), FALSE, TRUE)
      XS_tag[i] <- ifelse(all(grepl("\\d", (mcols(first(r_test))$XS))),
                          TRUE, FALSE)

    }


    # -- Setting the parameters for reading the BAM files

    bfl <- BamFileList(bam_files, asMates = TRUE, yieldSize = 100000)

    # If neither NH nor XS tags are available, the suboptimal alignment score
    # cannot be computed and the 3rd filter cannot be applied. Plus, unique
    # reads cannot be differentiated from multi-mapping reads
    if (!all(NH_tag) & !all(XS_tag)) {

      if (eval_unique_r == FALSE) {
        stop("Error in 'eval_unique_r = FALSE': Neither the NH tag nor the XS tag (from BWA) are provided in the BAM files. Unique reads cannot be differentiated from multi-mapping reads, therefore unique reads must also be filtered.")
      }

      warning("Neither the NH tag nor the XS tag (from BWA) are provided in the BAM files. At least one of the two tags is needed to obtain the suboptimal alignment score. The 3rd filter (based on the suboptimal alignment score) will not be applied.")
      param <- ScanBamParam(flag = sbflags,
                            tag = c("nM", "NM", "AS", "NH", "XS"))

    } else if (all(NH_tag) & !all(XS_tag)) {
      param <- ScanBamParam(flag = sbflags,
                            tag = c("nM", "NM", "AS", "NH", "XS"),
                            tagFilter = list(NH = c(2:10000)))

    } else if (!all(NH_tag) & all(XS_tag)) {
      param <- ScanBamParam(flag = sbflags,
                            tag = c("nM", "NM", "AS", "NH", "XS"),
                            tagFilter = list(XS = c(1:10000)))

    } else {
      stop("Error: an error occurred when reading the tags from the BAM files.")
    }

    table_counts <- matrix(nrow = length(ERV_ann), ncol = length(bfl))
    colnames(table_counts) <- gsub("\\.bam", "", names(bfl))


    for (i in 1:length(bfl)) {

      file <- bfl[[i]]
      open(file)

      while (length(r <- readGAlignmentPairs(file, use.names = TRUE,
                                             param = param,
                                             strandMode = strandMode))) {

        # Compute the suboptimal alignment score (equivalent to XS tag in BWA)
        # only if the score is not available in the BAM file and the NH tag is
        # available
        if (all(NH_tag) & !all(XS_tag)) {

          if (is.null(names(r))) {
            stop("The GAlignmentPairs object must have names")
          }

          r <- atena::getSuboptimalAlignScore(r, strandMode = strandMode)
        }

        r_first <- first(r)
        r_last <- last(r)
        rm(r)

        # The first filter is always applied
        cigar_first <- explodeCigarOpLengths(cigar(r_first), ops = c("H","S"))
        cigar_last <- explodeCigarOpLengths(cigar(r_last), ops = c("H","S"))

        SH_clipping_first <- lapply(cigar_first, function(x) sum(unlist(x)))
        SH_clipping_last <- lapply(cigar_last, function(x) sum(unlist(x)))


        # The second filter is always applied but when the NM tag is not
        # available the nM tag can be used instead.
        if (all(NM_tag)) {
          NM_filter_first <- (mcols(r_first)$NM / qwidth(r_first)) < 0.02
          NM_filter_last <- (mcols(r_last)$NM / qwidth(r_last)) < 0.02

        } else if (all(nM_tag)) {
          NM_filter_first <- (mcols(r_first)$nM / qwidth(r_first)) < 0.02
          NM_filter_last <- (mcols(r_last)$nM / qwidth(r_last)) < 0.02

        } else {
          stop("Error: Neither the NM tag nor the nM tag are available. At least one of them is needed to apply the 2nd filter.")
        }

        # The third filter (AS - XS >= 5) cannot be applied when neither the NH
        # tag nor the XS tag are available
        if (!all(NH_tag) & !all(XS_tag)) {
          AS_XS_filter_first <- rep(TRUE, length(r_first))
          AS_XS_filter_last <- rep(TRUE, length(r_last))
        } else {
          AS_XS_filter_first <- (mcols(r_first)$AS - mcols(r_first)$XS) >= 5
          AS_XS_filter_first[is.na(AS_XS_filter_first)] <- FALSE
          AS_XS_filter_last <- (mcols(r_last)$AS - mcols(r_last)$XS) >= 5
          AS_XS_filter_last[is.na(AS_XS_filter_last)] <- FALSE
        }

        # Filtering alignments considering all filters
        r_to_keep_f <- ((unlist(SH_clipping_first) / qwidth(r_first)) < 0.02) &
          NM_filter_first & AS_XS_filter_first

        r_to_keep_l <- ((unlist(SH_clipping_last) / qwidth(r_last)) < 0.02) &
          NM_filter_last & AS_XS_filter_last


        # Selecting the filtered fragments, considering the two fragments of a
        # paired-end read as independent and Joining the fragments into a
        # unique object
        r_total <- c(r_first[r_to_keep_f], r_last[r_to_keep_l])
        rm(r_first, r_last)

        # Computing read counts for each element
        overlap <- summarizeOverlaps(features = ERV_ann, reads = r_total,
                                     ignore.strand = ignore_strand,
                                     mode = Union, inter.feature = FALSE)
        rm(r_total)

        if (all(is.na(table_counts[,i]))) {
          table_counts[,i] <- assay(overlap)
        } else {
          table_counts[,i] <- table_counts[,i] + assay(overlap)
        }

      }

      close(file)

      if (any(all(NH_tag) & all(XS_tag))) {
        ## If the NH and XS are not provided, in the previous step both unique
        ## reads and multimapping reads have already been filtered
        ## If not, if the user wants to filter unique reads, they will be
        ## filtered (only the first 2 filters). Instead, if the user does not
        ## want to filter unique reads, they will not be filtered, but they
        ## will be considered for computing HERV raw counts.

        # Now the alignments corresponding to unique reads are read
        sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                               isSupplementaryAlignment = FALSE,
                               isProperPair = TRUE,
                               isDuplicate = FALSE,
                               isNotPassingQualityControls = FALSE)


        if (!all(NH_tag) & all(XS_tag)) {

          # If the NH is not available but the XS is, only unique reads can be
          # read by reading those in which XS = 0.
          param <- ScanBamParam(flag = sbflags,
                                tag = c("nM", "NM", "AS", "NH", "XS"),
                                tagFilter = list(XS = 0))
        } else {
          param <- ScanBamParam(flag = sbflags,
                                tag = c("nM","NM","AS","NH","XS"),
                                tagFilter = list(NH = 1)) # Read only unique reads
        }

        open(file)

        while (length(r <- readGAlignmentPairs(file, param = param,
                                               strandMode = strandMode))) {
          r_first <- first(r)
          r_last <- last(r)
          rm(r)

          # If the user wants the unique reads to pass the filters in ERVmap,
          # run ERVmap (except for the "AS - XS > 5" filter) on the unique reads
          if (eval_unique_r) {

            cigar_first <- explodeCigarOpLengths(cigar(r_first),
                                                 ops =c("H","S"))
            cigar_last <- explodeCigarOpLengths(cigar(r_first),
                                                ops = c("H", "S"))

            SH_clipping_first <- lapply(cigar_first, function(x) sum(unlist(x)))
            SH_clipping_last <- lapply(cigar_last, function(x) sum(unlist(x)))

            if (all(NM_tag)) {
              NM_filter_first <- (mcols(r_first)$NM / qwidth(r_first)) < 0.02
              NM_filter_last <- (mcols(r_last)$NM / qwidth(r_last)) < 0.02

            } else if (all(nM_tag)) {
              NM_filter_first <- (mcols(r_first)$nM / qwidth(r_first)) < 0.02
              NM_filter_last <- (mcols(r_last)$nM / qwidth(r_last)) < 0.02

            } else {
              stop("Error: Neither the NM tag nor the nM tag are available. At least one of them is needed to apply the 2nd filter.")
            }

            # Filtering alignments considering the two filters
            r_to_keep_f <- ((unlist(SH_clipping_first) / qwidth(r_first)) < 0.02) &
              NM_filter_first

            r_to_keep_l <- ((unlist(SH_clipping_last) / qwidth(r_last)) < 0.02) &
              NM_filter_last

            r_unique <- c(r_first[r_to_keep_f], r_last[r_to_keep_l])

            # If the user does not want the unique reads to be filtered, but
            # instead wants to include all unique alignments in the analysis
            # (discarding unmapped, unpaired, duplicated and not passing
            # quality control reads), then:
          } else {
              r_unique <- c(r_first, r_last)
          }

          # Computing read counts for each element
          overlap <- summarizeOverlaps(features = ERV_ann, reads = r_unique,
                                       ignore.strand = ignore_strand,
                                       mode = Union, inter.feature = FALSE)
          rm(r_first, r_last, r_unique)

          if (all(is.na(table_counts[,i]))) {
            table_counts[,i] <- assay(overlap)
          } else {
            table_counts[,i] <- table_counts[,i] + assay(overlap)
          }

        }

        close(file)
        gc()

      }

      sample_name <- gsub("\\.bam", "", names(bfl)[i])
      colnames(table_counts)[i] <- sample_name
      message(paste("Sample", sample_name, "done!", sep = " "))

    }
  }
  df_coord <- data.frame(chr =  seqnames(ERV_ann),
                         start = start(ranges(ERV_ann)),
                         end = end(ranges(ERV_ann)),
                         strand = strand(ERV_ann))

  table_counts_coord <- cbind(df_coord, table_counts)
  rownames(table_counts_coord) <- rowData(overlap)$name

  return(SummarizedExperiment::makeSummarizedExperimentFromDataFrame(table_counts_coord))

}
