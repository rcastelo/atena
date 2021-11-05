test_TEtranscripts <- function() {
    bamfiles <- list.files(system.file("extdata", package="atena"),
                                pattern="*.bam", full.names=TRUE)
    TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                                package="atena"))
    ttpar <- TEtranscriptsParam(bamfiles, teFeatures = TE_annot,
                                    singleEnd = TRUE, ignoreStrand=TRUE,
                                    aggregateby = c("repName"))
    ttSE <- qtex(ttpar)
    
    checkEqualsNumeric(dim(ttSE), c(1, 2))
    checkEqualsNumeric(head(sort(assay(ttSE), decreasing=TRUE)), 
                            c(149, 122))
}


test_tt_input1 <- function() {
    
    TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                            package="atena"))
    checkException(TEtranscriptsParam(teFeatures=TE_annot, singleEnd = TRUE),
                    "An error prompts when no input BAM file is specified",
                    silent=TRUE)
}


test_tt_input2 <- function() {
    
    bamfiles <- list.files(system.file("extdata", package="atena"),
                            pattern="*.bam", full.names=TRUE)
    checkException(TEtranscriptsParam(bfl=bamfiles, singleEnd = TRUE),
                    "An error prompts when TE annotations are not specified",
                    silent=TRUE)
}


test_exon_summarization <- function() {
    
    bamfiles <- list.files(system.file("extdata", package="atena"),
                            pattern="*.bam", full.names=TRUE)
    TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                            package="atena"))
    
    # Creating an example of gene annotations
    annot_gen <- GRanges(seqnames = rep("2L",10),
                        ranges = IRanges(
                            start = c(1,20,45,80,110,130,150,170,200,220),
                            width = c(10,20,35,10,5,15,10,25,5,20)),
                        strand = "*", 
                        type = rep("exon",10))
    # Setting gene ids
    names(annot_gen) <- paste0("gene",c(rep(1,3),rep(2,4),rep(3,3)))
    annot_gen
    ttpar_gen <- TEtranscriptsParam(bamfiles, teFeatures = TE_annot,
                                    geneFeatures = annot_gen, singleEnd = TRUE,
                                    ignoreStrand=TRUE)
    
    checkTrue(is(features(ttpar_gen), "GRangesList"),
                "Exons are summarized at the gene level in a GRangesList
                object")
}
