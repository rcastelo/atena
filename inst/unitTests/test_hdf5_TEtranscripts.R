test_hdf5_TEtranscripts <- function() {
    bamfiles <- list.files(system.file("extdata", package="atena"),
                                pattern="*.bam", full.names=TRUE)
    TE_annot <- readRDS(file = system.file("extdata", "Top28TEs.rds",
                                                package="atena"))
    ttpar <- TEtranscriptsParam(bamfiles, teFeatures = TE_annot,
                                    singleEnd = TRUE, ignoreStrand=TRUE,
                                    aggregateby = c("repName"))
    ttSE1 <- qtex(ttpar)
	ttSE2 <- qtex(ttpar, on.disk=TRUE)
    
    checkEqualsNumeric(dim(ttSE1), c(1, 2))
    checkEqualsNumeric(head(sort(assay(ttSE1), decreasing=TRUE)),
                            c(149, 122))
	
    checkEqualsNumeric(dim(ttSE1), c(1, 2))
    checkEqualsNumeric(head(sort(assay(ttSE1), decreasing=TRUE)),
                            c(149, 122))
}
