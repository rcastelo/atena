test_ERVmap <- function() {
    bamfiles <- list.files(system.file("extdata", package="atena"),
                                pattern="*.bam", full.names=TRUE)
    TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                            package="atena"))
    empar <- ERVmapParam(bamfiles, teFeatures = TE_annot, singleEnd = TRUE,
                            ignoreStrand = TRUE, suboptimalAlignmentCutoff=NA)
    emSE <- qtex(empar)
    
    checkEqualsNumeric(dim(emSE), c(28, 2))
    checkEqualsNumeric(head(sort(assay(emSE), decreasing=TRUE)), 
                            c(4, 3, 3, 3, 2, 2))
}

test_em_input1 <- function() {
    
    TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                            package="atena"))
    checkException(ERVmapParam(teFeatures=TE_annot, singleEnd = TRUE),
                    "An error prompts when no input BAM file is specified",
                    silent=TRUE)
}


test_em_input2 <- function() {
    
    bamfiles <- list.files(system.file("extdata", package="atena"),
                            pattern="*.bam", full.names=TRUE)
    checkException(ERVmapParam(bfl=bamfiles, singleEnd = TRUE),
                    "An error prompts when TE annotations are not specified",
                    silent=TRUE)
}



test_combine_seqlevels_annot <- function() {
    
    bamfiles <- list.files(system.file("extdata", package="atena"),
                            pattern="*.bam", full.names=TRUE)
    TE_annot <- GRanges(seqnames = rep("3R",8),
                        ranges = IRanges(
                            start = c(30,42,84,120,134,159,175,250),
                            width = c(35,10,5,15,10,25,5,20)),
                        strand = "*")
    names(TE_annot) <- paste("TE", seq(1,8), sep = "")
    gene_annot <- GRanges(seqnames = rep("2L",10),
                        ranges = IRanges(
                            start = c(1,20,45,80,110,130,150,170,200,220),
                            width = c(10,20,35,10,5,15,10,25,5,20)),
                        strand = "*", 
                        type = rep("exon",10))
    names(gene_annot) <- paste("gene", seq(1,10), sep = "")
    empar <- ERVmapParam(bamfiles, teFeatures = TE_annot, 
                        geneFeatures = gene_annot, singleEnd = TRUE,
                        ignoreStrand = TRUE, suboptimalAlignmentCutoff=NA)
    
    checkTrue(seqlevels(TE_annot) %in% seqlevels(features(empar)))
    checkTrue(seqlevels(gene_annot) %in% seqlevels(features(empar)))
}




