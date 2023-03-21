test_Telescope <- function() {
    bamfiles <- list.files(system.file("extdata", package="atena"),
                                pattern="*.bam", full.names=TRUE)
    TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                            package="atena"))
    gene_annot <- readRDS(file = system.file("extdata", "Top50genes.rds",
                                                package="atena"))
    tspar <- TelescopeParam(bfl=bamfiles, teFeatures=TE_annot,
                                geneFeatures = gene_annot,
                                singleEnd = TRUE, ignoreStrand=TRUE)
    tsSE <- qtex(tspar)
    
    checkEqualsNumeric(dim(tsSE), c(76, 2))
    checkEqualsNumeric(head(sort(assay(tsSE), decreasing=TRUE)), 
                            as.integer(c(149, 144, 6, 1, 0, 0)))  # c(150, 150, 0, 0, 0, 0)

}


test_ts_input1 <- function() {
    
    TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                            package="atena"))
    checkException(TelescopeParam(teFeatures=TE_annot, singleEnd = TRUE),
                    "An error prompts when no input BAM file is specified",
                    silent=TRUE)
}


test_ts_input2 <- function() {
    
    bamfiles <- list.files(system.file("extdata", package="atena"),
                           pattern="*.bam", full.names=TRUE)
    checkException(TelescopeParam(bfl=bamfiles, singleEnd = TRUE),
                    "An error prompts when TE annotations are not specified",
                    silent=TRUE)
}





