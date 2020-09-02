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
#'   class for further detail.
#'
#' @param ann_file (Default NULL) A single string containing the path and name
#'   of the endogenous retrovirus (ERVs) or TE annotation file (including the
#'   extension of the file). The allowed file formats are BED and GTF / GFF
#'   version 2. When \code{ann_file} is not specified (\code{ann_file} == NULL),
#'   the human endogenous retroviruses annotations published by
#'   \href{https://doi.org/10.1073/pnas.1814589115}{Tokuyama et al. (2018)}
#'   will be used instead.
#'
#' @param eval_unique_r (Default TRUE) Logical value indicating whether to apply
#' the alignment filters to unique reads (TRUE) or not (FALSE). These filters,
#' which are always applied to multi-mapping reads, are optional for unique
#' reads. If TRUE, the unique reads not passing one or more filters from the
#' ERVmap pipeline, except for the "AS - XS >= 5" filter, will be discarded to
#' compute TEs expression.
#'
#' @param ignore_strand (Default FALSE) A logical which defines if the strand
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
#'
#' @return A `SummarizedExperiment` object representing containing the read
#' counts for each sample and element in the annotations file. Each sample is 
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
#' \dontrun{path_bam_files_sorted <- "Documents/RNAseq/BAMfiles_sorted"
#'
#' table_counts <- count_TE(path_bam_files_sorted, strandMode = 2,
#'                         eval_unique_r = TRUE,
#'                         ignore_strand = FALSE)}
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



count_TE <- function(path_bam_files_sorted, strandMode = 1,
                     ann_file = NULL, eval_unique_r = TRUE,
                     ignore_strand = FALSE) {

  ERVmap_ann <- r_XS <- reads_to_keep <- NULL

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

    } else {
      if (ERV_format == "gtf" | ERV_format == "gff") {
        ERV_ann <- rtracklayer::readGFF(ann_file)
        ERV_ann_granges <- GenomicRanges::makeGRangesFromDataFrame(ERV_ann)
        mcols(ERV_ann_granges) <- ERV_ann$gene_id
        colnames(mcols(ERV_ann_granges))[1] <- "name"
        ERV_ann <- ERV_ann_granges
      }
    }

    # Changing "seqnames" to match the ones from the reads
    if (!all(grepl(x = levels(seqnames(ERV_ann)),
                   pattern = "chr", fixed = TRUE))) {
      GenomeInfoDb::seqlevels(ERV_ann) <- levels(as.factor(paste("chr",
                                                                 seqnames(ERV_ann),
                                                   sep = '')))
    }
  }


  # -- Reading the BAM files
  bam_files <- list.files(path_bam_files_sorted, pattern = "\\.bam$",
                          full.names = TRUE)
  bfl <- BamFileList(bam_files, asMates = TRUE, yieldSize = 10000000)

  sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                         isProperPair = TRUE,
                         isDuplicate = FALSE,
                         isNotPassingQualityControls = FALSE)

  param <- ScanBamParam(flag = sbflags,
                        tag = c("nM", "NM", "AS", "NH", "XS"),
                        what = scanBamWhat(),
                        tagFilter = list(NH = c(2:10000)))


  table_counts <- matrix(nrow = length(ERV_ann), ncol = length(bfl))


  for (i in 1:length(bfl)) {

    file <- bfl[[i]]
    open(file)

    r_multi_filtered <- NULL

    while (length(r <- readGAlignmentPairs(file, use.names=TRUE,
                                           param=param,
                                           strandMode=strandMode))) {

      if (is.null(names(r))) {
        stop("The GAlignmentPairs object must have names")
      }

      r <- atena::getSuboptimalAlignScore(r, strandMode = strandMode)

      # Now that the XS tag has been obtained for the multi-mapping reads, the
      # ERVmap pipeline will be applied to these alignments

      cigar_first <- explodeCigarOpLengths(cigar(first(r_XS)), ops = c("H","S"))
      cigar_last <- explodeCigarOpLengths(cigar(last(r_XS)), ops = c("H","S"))

      SH_clipping_first <- lapply(cigar_first, function(x) sum(unlist(x)))
      SH_clipping_last <- lapply(cigar_last, function(x) sum(unlist(x)))

      # AS - XS >= 5 filter
      AS_XS_filter <- (mcols(first(r_XS))$AS - mcols(first(r_XS))$XS) >= 5 |
        (mcols(last(r_XS))$AS - mcols(last(r_XS))$XS) >= 5
      AS_XS_filter[is.na(AS_XS_filter)] <- FALSE

      # Filtering alignments considering all filters
      r_to_keep <- ((unlist(SH_clipping_first) / qwidth(first(r_XS))) < 0.02 |
                      (unlist(SH_clipping_last) / qwidth(last(r_XS))) < 0.02 ) &
        ((mcols(first(r_XS))$NM / qwidth(first(r_XS))) < 0.02 |
           (mcols(last(r_XS))$NM / qwidth(last(r_XS))) < 0.02)  &
        AS_XS_filter

      if (is.null(r_multi_filtered)) {
        r_multi_filtered <- r_XS[r_to_keep]
      }
      else {
        r_multi_filtered <- c(r_multi_filtered, r_XS[r_to_keep])
      }

    }


    # Now the alignments corresponding to unique reads are read
    sbflags <- scanBamFlag(isUnmappedQuery = FALSE,
                           isSupplementaryAlignment = FALSE,
                           isProperPair = TRUE,
                           isDuplicate = FALSE,
                           isNotPassingQualityControls = FALSE)

    param <- ScanBamParam(flag = sbflags,
                          tag = c("nM","NM","AS","NH","XS"),
                          tagFilter = list(NH = 1)) # Read only unique reads

    r_unique <- NULL

    while (length(r <- readGAlignmentPairs(file, use.names = TRUE,
                                           param = param,
                                           strandMode = strandMode))) {

      # If the user wants the unique reads to pass the filters in ERVmap, run
      # ERVmap (except for the "AS - XS > 5" filter) on the unique reads
      if (eval_unique_r) {

        cigar_first <- explodeCigarOpLengths(cigar(first(r_XS)), ops =c("H","S"))
        cigar_last <- explodeCigarOpLengths(cigar(last(r_XS)), ops = c("H", "S"))

        SH_clipping_first <- lapply(cigar_first, function(x) sum(unlist(x)))
        SH_clipping_last <- lapply(cigar_last, function(x) sum(unlist(x)))

        r_to_keep <- ((unlist(SH_clipping_first) / qwidth(first(r_XS))) < 0.02 |
                        (unlist(SH_clipping_last)/ qwidth(last(r_XS))) < 0.02) &
          ((mcols(first(r_XS))$NM / qwidth(first(r_XS))) < 0.02 |
             (mcols(last(r_XS))$NM / qwidth(last(r_XS))) < 0.02)

        if (is.null(r_unique)) {
          r_unique <- r[reads_to_keep]
        }
        else {
          r_unique <- c(r_unique, r[reads_to_keep])
        }
      }

      # If the user does not want to select only the unique reads that pass the
      # ERVmap filters, but instead wants to include all unique alignments in
      # the analysis (discarding unmapped, unpaired, duplicated and not passing
      # quality control reads), then:

      else {
        if (is.null(r_unique)) {
          r_unique <- r
        }
        else {
          r_unique <- c(r_unique,r)
        }

      }

    }

    close(file)

    # Joining the alignments corresponding to multi-mapping and unique reads
    r_total <- c(r_multi_filtered, r_unique)

    overlap <- summarizeOverlaps(features = ERV_ann, reads = r_total,
                                 ignore.strand = ignore_strand,
                                 strandMode=strandMode, mode = Union,
                                 inter.feature = FALSE, singleEnd = FALSE)

    table_counts[,i] <- assay(overlap)

    sample_name <- gsub("\\.bam", "",
                        grep("\\.bam", strsplit(path(file), "/")[[1]],
                             value = TRUE))
    colnames(table_counts)[i] <- sample_name

    print(paste("Sample", sample_name, "done!", sep = " "))

  }
  
  
  df_coord <- data.frame(chr =  seqnames(ERV_ann), 
                          start = start(ranges(ERV_ann)), 
                          end = end(ranges(ERV_ann)),
                          strand = strand(ERV_ann))
  
  table_counts_coord <- cbind(df_coord, table_counts)
  rownames(table_counts_coord) <- rowData(overlap)$name
  
  return(SummarizedExperiment::makeSummarizedExperimentFromDataFrame(table_counts_coord))

}
