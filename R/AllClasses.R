#' Atena parameter class
#'
#' This is a virtual class from which other classes are derived
#' for storing parameters provided to quantification methods of
#' transposable elements from RNA-seq data.
#'
#' @slot bfl A \linkS4class{BamFileList} object.
#'
#' @slot features A \linkS4class{GRanges} object.
#'
#' @slot aggregateby Character vector with column names in the annotation
#'       to be used to aggregate quantifications.
#'       
#' @slot ovMode Character vector indicating the overlapping mode. Available
#'       options are: "ovUnion" (default) and "ovIntersectionStrict",
#'       which implement the corresponding methods from HTSeq
#'       (\url{https://htseq.readthedocs.io/en/release_0.11.1/count.html}).
#'       In the TEtranscripts, ERVmap and Telescope methods ambiguous
#'       alignments (alignments overlapping > 1 feature) are addressed
#'       differently depending on the method. In the atena method, those
#'       overlaps are not counted.
#'    
#'
#' @importClassesFrom Rsamtools BamFileList
#' @importClassesFrom GenomicRanges GenomicRanges_OR_GenomicRangesList
#'
#' @seealso
#' \code{\link{ERVmapParam-class}}
#' \code{\link{TelescopeParam-class}}
#' \code{\link{TEtranscriptsParam-class}}
#' \code{\link{atenaParam-class}}
#'
#' @examples
#' bamfiles <- list.files(system.file("extdata", package="atena"),
#'                        pattern="*.bam", full.names=TRUE)
#' TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
#'                     package="atena"))
#' ttpar <- TEtranscriptsParam(bamfiles, teFeatures=TE_annot, singleEnd=TRUE,
#'                             ignoreStrand=TRUE, aggregateby = c("repName"))
#' path(ttpar)
#'
#' @name QuantifyParam-class
#' @rdname QuantifyParam-class
#' @exportClass QuantifyParam
setClass("QuantifyParam",
        representation(bfl="BamFileList",
                        features="GenomicRanges_OR_GenomicRangesList",
                        aggregateby="character", ovMode="character"))

#' @param object A \linkS4class{QuantifyParam} object.
#'
#' @importFrom BiocGenerics path
#'
#' @return \code{path()}: Filesystem paths to the BAM files in the input
#' parameter object.
#'
#' @export
#' @aliases path,QuantifyParam-method
#' @rdname QuantifyParam-class
setMethod("path", "QuantifyParam",
        function(object) {
            path(object@bfl)
        })


#' @param object A \linkS4class{QuantifyParam} object.
#'
#' @return \code{features()}: The \code{GenomicRanges} or
#' \code{GenomicRangesList} object with the features in the input parameter
#' object.
#'
#' @export
#' @aliases features
#' @aliases features,QuantifyParam-method
#' @rdname QuantifyParam-class
#' @include AllGenerics.R
setMethod("features", "QuantifyParam",
        function(object) {
            object@features
        })

#' ERVmap parameter class
#'
#' This is a class for storing parameters provided to the ERVmap algorithm.
#' It is a subclass of the 'QuantifyParam-class'.
#'
#' @slot readMapper The name of the software used to align reads, obtained from
#' the BAM file header.
#'
#' @slot singleEnd (Default FALSE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @slot strandMode (Default 1) Numeric vector which can take values 0, 1 or 2.
#'   The strand mode is a per-object switch on
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   objects that controls the behavior of the strand getter. See
#'   \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#'   class for further detail. If \code{singleEnd = TRUE}, then
#'   \code{strandMode} #'   is ignored.
#'
#' @slot ignoreStrand (Default TRUE) A logical which defines if the strand
#' should be taken into consideration when computing the overlap between reads
#' and TEs in the annotations. When \code{ignore_strand = FALSE}, only those
#' reads which overlap the TE and are on the same strand are counted. On the
#' contrary, when \code{ignore_strand = TRUE}, any read overlapping an element
#' in \code{teFeatures} is counted regardless of the strand.
#'
#' @slot fragments (Default not \code{singleEnd}) A logical; applied to
#' paired-end data only. When \code{fragments=TRUE}, the read-counting
#' method in the original ERVmap algorithm is applied: each mate of a
#' paired-end read is counted (including ambiguous and not properly paired 
#' reads). When
#' \code{fragments=FALSE}, if the two mates of a paired-end read map to the
#' same element, they are counted as a single hit and singletons, reads with
#' unmapped pairs and other ambiguous or not properly paired fragments are 
#' not counted (see "Pairing criteria" in 
#' \code{\link[GenomicAlignments]{readGAlignments}()}).
#'
#' @slot maxMismatchRate (Default 0.02) Numeric value storing the maximum
#' mismatch rate employed by the ERVmap algorithm to discard aligned reads
#' whose rate of sum of hard and soft clipping, or of the edit distance over
#' the genome reference, to the length of the read is above this threshold.
#'
#' @slot suboptimalAlignmentTag (Default "auto") Character string storing the
#' tag name in the BAM files that stores the suboptimal alignment score used in
#' the third filter of ERVmap; see Tokuyama et al. (2018). The default,
#' \code{suboptimalAlignmentTag="auto"}, assumes that either the BAM files were
#' generated by BWA and include a tag called \code{XS} that stores the
#' suboptimal alignment score or, if the \code{XS} tag is not available, then
#' it uses the available secondary alignments to implement an analogous
#' approach to that of the third ERVmap filter. When
#' \code{suboptimalAlignmentTag="none"}, it also performs the latter approach
#' even when the tag \code{XS} is available.
#' When this parameter is different from \code{"auto"} and \code{"none"}, a tag
#' with the given name is used to extract the suboptimal alignment score.
#' The absence of that tag will prompt an error.
#'
#' @slot suboptimalAlignmentCutoff (Default 5) Numeric value storing the cutoff
#' above which the difference between the alignment score and the suboptimal
#' alignment score is considered sufficiently large to retain the alignment.
#' When this value is set to \code{NA}, then the filtering step based on
#' suboptimal alignment scores is skipped.
#'
#' @slot geneCountMode (Default "all") Character string indicating if the
#' ERVmap read filters applied to quantify TEs expression should also be
#' applied when quantifying gene expression ("ervmap") or not ("all"), in which
#' case all primary alignments mapping to genes are counted.
#'
#' @references
#' Tokuyama M et al. ERVmap analysis reveals genome-wide transcription of human
#' endogenous retroviruses. PNAS. 2018;115(50):12565-12572. DOI:
#' \url{https://doi.org/10.1073/pnas.1814589115}
#'
#' @name ERVmapParam-class
#' @rdname ERVmapParam-class
#' @exportClass ERVmapParam
setClass("ERVmapParam", contains="QuantifyParam",
        representation(singleEnd="logical",
                        ignoreStrand="logical",
                        strandMode="integer",
                        fragments="logical",
                        maxMismatchRate="numeric",
                        readMapper="character",
                        suboptimalAlignmentTag="character",
                        suboptimalAlignmentCutoff="numeric",
                        geneCountMode="character"))

#' Telescope parameter class
#'
#' This is a class for storing parameters
#' provided to the Telescope algorithm.
#'
#' @slot singleEnd (Default TRUE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @slot strandMode (Default 1) Numeric vector which can take values 0, 1 or 2.
#' The strand mode is a per-object switch on
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' objects that controls the behavior of the strand getter. See
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' class for further detail. If \code{singleEnd = TRUE}, then \code{strandMode}
#' is ignored.
#'
#' @slot ignoreStrand (Default FALSE) A logical which defines if the strand
#' should be taken into consideration when computing the overlap between reads
#' and annotated features. When \code{ignoreStrand = FALSE}, an aligned read
#' is considered to be overlapping an annotated feature as long as they
#' have a non-empty intersecting genomic range on the same strand, while when
#' \code{ignoreStrand = TRUE} the strand is not considered.
#'
#' @slot fragments (Default FALSE) A logical; applied to paired-end data only.
#' When \code{fragments=FALSE}, the read-counting method only counts
#' ‘mated pairs’ from opposite strands (non-ambiguous properly paired reads), 
#' while when \code{fragments=TRUE} same-strand pairs, singletons, reads with 
#' unmapped pairs and other ambiguous or not properly paired fragments
#' are also counted (see "Pairing criteria" in 
#' \code{\link[GenomicAlignments]{readGAlignments}()}). \code{fragments=TRUE} 
#' is equivalent to the original Telescope algorithm. For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}()}.
#' 
#' @slot minOverlFract (Default 0.2) A numeric scalar. \code{minOverlFract}
#' is multiplied by the read length and the resulting value is used to
#' discard alignments for which the overlapping length (number of base
#' pairs the alignment and the feature overlap) is lower. When no minimum 
#' overlap is required, set \code{minOverlFract = 0}.
#'
#' @slot pi_prior (Default 0) A positive integer scalar indicating the prior
#' on pi. This is equivalent to adding n unique reads.
#'
#' @slot theta_prior (Default 0) A positive integer scalar storing the prior
#' on Q. Equivalent to adding n non-unique reads.
#'
#' @slot em_epsilon (Default 1e-7) A numeric scalar indicating the EM
#' Algorithm Epsilon cutoff.
#'
#' @slot maxIter A positive integer scalar storing the maximum number of
#' iterations of the EM SQUAREM algorithm (Du and Varadhan, 2020). Default
#' is 100 and this value is passed to the \code{maxiter} parameter of the
#' \code{\link[SQUAREM]{squarem}()} function.
#' 
#' @slot reassign_mode (Default 'exclude') Character vector indicating
#' reassignment mode after EM step. 
#' Available methods are 'exclude' (reads with more than one best
#' assignment are excluded from the final counts), 'choose' (when reads have
#' more than one best assignment, one of them is randomly chosen), 'average'
#' (the read count is divided evenly among the best assignments) and 'conf'
#' (only assignments that exceed a certain threshold -defined by 
#' \code{conf_prob} parameter- are accepted, then the read count is
#' proportionally divided among the assignments above \code{conf_prob}).
#' 
#' @slot conf_prob (Default 0.9) Minimum probability for high confidence
#' assignment.
#' 
#' 
#' @references
#' Bendall et al. Telescope: characterization of the retrotranscriptome by
#' accurate estimation of transposable element expression.
#' PLOS Comp. Biol. 2019;15(9):e1006453. DOI:
#' \url{https://doi.org/10.1371/journal.pcbi.1006453}
#'
#' @name TelescopeParam-class
#' @rdname TelescopeParam-class
#' @exportClass TelescopeParam
setClass("TelescopeParam", contains="QuantifyParam",
        representation(singleEnd="logical",
                        strandMode="integer",
                        ignoreStrand="logical",
                        fragments="logical",
                        minOverlFract="numeric",
                        pi_prior="integer",
                        theta_prior="integer",
                        em_epsilon="numeric",
                        maxIter="integer",
                        reassign_mode="character",
                        conf_prob="numeric"))

#' TEtranscripts parameter class
#'
#' This is a class for storing parameters provided to the TEtranscripts
#' algorithm. It is a subclass of the 'QuantifyParam-class'.
#'
#' @slot singleEnd (Default FALSE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @slot ignoreStrand (Default FALSE) A logical which defines if the strand
#' should be taken into consideration when computing the overlap between reads
#' and annotated features. When \code{ignoreStrand = FALSE}, an aligned read
#' will be considered to be overlapping an annotated feature as long as they
#' have a non-empty intersecting genomic ranges on the same strand, while when
#' \code{ignoreStrand = TRUE} the strand will not be considered.
#'
#' @slot strandMode (Default 1) Numeric vector which can take values 0, 1 or 2.
#' The strand mode is a per-object switch on
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' objects that controls the behavior of the strand getter. See
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' class for further detail. If \code{singleEnd = TRUE}, then use either
#' \code{strandMode = NULL} or do not specify the \code{strandMode} parameter.
#'
#' @slot fragments (Default TRUE) A logical; applied to paired-end data only.
#' In both cases (\code{fragments=FALSE} and \code{fragments=TRUE}), the
#' read-counting method discards not properly paired reads. Moreover,
#' when \code{fragments=FALSE}, only non-ambiguous properly paired reads are
#' counted. When \code{fragments=TRUE}, ambiguous reads are also counted 
#' (see "Pairing criteria" in \code{\link[GenomicAlignments]{readGAlignments}()}). 
#' \code{fragments=TRUE} is equivalent to
#' the behavior of the TEtranscripts algorithm. For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}()}.
#'
#' @slot tolerance A positive numeric scalar storing the minimum tolerance
#' above which the SQUAREM algorithm (Du and Varadhan, 2020) keeps iterating.
#' Default is \code{1e-4} and this value is passed to the \code{tol} parameter
#' of the \code{\link[SQUAREM]{squarem}()} function.
#'
#' @slot maxIter A positive integer scalar storing the maximum number of
#' iterations of the SQUAREM algorithm (Du and Varadhan, 2020). Default
#' is 100 and this value is passed to the \code{maxiter} parameter of the
#' \code{\link[SQUAREM]{squarem}()} function.
#'
#' @references
#' Jin Y et al. TEtranscripts: a package for including transposable elements
#' in differential expression analysis of RNA-seq datasets.
#' Bioinformatics. 2015;31(22):3593-3599. DOI:
#' \url{https://doi.org/10.1093/bioinformatics/btv422}
#'
#' @name TEtranscriptsParam-class
#' @rdname TEtranscriptsParam-class
#' @exportClass TEtranscriptsParam
setClass("TEtranscriptsParam", contains="QuantifyParam",
        representation(singleEnd="logical",
                        ignoreStrand="logical",
                        strandMode="integer",
                        fragments="logical",
                        tolerance="numeric",
                        maxIter="integer"))

#' atena parameter class
#'
#' This is a class for storing parameters to quantify TE (and gene) expression
#' using the atena method. It is a subclass of the 'QuantifyParam-class'.
#'
#' @slot singleEnd (Default TRUE) Logical value indicating if reads are single
#' (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @slot strandMode (Default 1) Numeric vector which can take values 0, 1 or 2.
#' The strand mode is a per-object switch on
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' objects that controls the behavior of the strand getter. See
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' class for further detail. If \code{singleEnd = TRUE}, then \code{strandMode}
#' is ignored.
#'
#' @slot ignoreStrand (Default FALSE) A logical which defines if the strand
#' should be taken into consideration when computing the overlap between reads
#' and annotated features. When \code{ignoreStrand = FALSE}, an aligned read
#' is considered to be overlapping an annotated feature as long as they
#' have a non-empty intersecting genomic range on the same strand, while when
#' \code{ignoreStrand = TRUE} the strand is not considered.
#'
#' @slot fragments (Default TRUE) A logical; applied to paired-end data only.
#' When \code{fragments=FALSE}, the read-counting method only counts
#' ‘mated pairs’ from opposite strands (non-ambiguous properly paired reads), 
#' while when \code{fragments=TRUE} same-strand pairs, singletons, reads with 
#' unmapped pairs and other ambiguous or not properly paired fragments
#' are also counted (see "Pairing criteria" in 
#' \code{\link[GenomicAlignments]{readGAlignments}()}). 
#' For further details see
#' \code{\link[GenomicAlignments]{summarizeOverlaps}()}.
#'
#' @slot pi_prior (Default 0) A positive numeric object indicating the prior
#' on pi. The same prior can be specified for all features setting
#' \code{pi_prior} as a scalar, or each feature can have a specific prior by
#' setting \code{pi_prior} as a vector with \code{names()} corresponding to
#' all feature names. Setting a pi prior is equivalent to adding n unique
#' reads.
#'
#' @slot theta_prior (Default 0) A positive numeric object indicating the
#' prior on Q. The same prior can be specified for all features setting
#' \code{theta_prior} as a scalar, or each feature can have a specific prior by
#' setting \code{theta_prior} as a vector with \code{names()} corresponding to
#' all feature names. Equivalent to adding n non-unique reads.
#'
#' @slot em_epsilon (Default 1e-7) A numeric scalar indicating the EM
#' Algorithm Epsilon cutoff.
#'
#' @slot maxIter A positive integer scalar storing the maximum number of
#' iterations of the EM SQUAREM algorithm (Du and Varadhan, 2020). Default
#' is 100 and this value is passed to the \code{maxiter} parameter of the
#' \code{\link[SQUAREM]{squarem}()} function.
#' 
#' @slot reassign_mode (Default 'exclude') Character vector indicating
#' reassignment mode after EM step. 
#' Available methods are 'exclude' (reads with more than one best
#' assignment are excluded from the final counts), 'choose' (when reads have
#' more than one best assignment, one of them is randomly chosen), 'average'
#' (the read count is divided evenly among the best assignments) and 'conf'
#' (only assignments that exceed a certain threshold -defined by 
#' \code{conf_prob} parameter- are accepted, then the read count is
#' proportionally divided among the assignments above \code{conf_prob}).
#' 
#' @slot conf_prob (Default 0.9) Minimum probability for high confidence
#' assignment.
#' 
#'
#' @name atenaParam-class
#' @rdname atenaParam-class
#' @exportClass atenaParam
setClass("atenaParam", contains="QuantifyParam",
         representation(singleEnd="logical",
                        strandMode="integer",
                        ignoreStrand="logical",
                        fragments="logical",
                        pi_prior="numeric",
                        theta_prior="numeric",
                        em_epsilon="numeric",
                        maxIter="integer",
                        reassign_mode="character",
                        conf_prob="numeric"))
