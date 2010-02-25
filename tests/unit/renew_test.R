sp <- SolexaPath(system.file("extdata", package="ShortRead"))
ap <- analysisPath(sp)
filt <- chromosomeFilter("chr[[:digit:]+].fa")
aln <- readAligned(ap, "s_2_export.txt", "SolexaExport",
                   filter=filt)

test_renewable0 <- function()
{
    cls <- renewable()
    for (cl in cls) {
        def <- getClass(cl, where=getNamespace("ShortRead"))
        checkTrue(validObject(df))
    }
}

test_renewable_non_virtual<- function()
{
    cls <- renewable()
    for (cl in cls) {
        if (!getClass(cl)@virtual)
            checkIdentical(getSlots(cl), renewable(cl)[[1]])
    }
}

test_renew <- function()
{
    checkIdentical(aln, renew(aln))

    labels <- sub("\\.fa", "", levels(chromosome(aln)))
    updt <- factor(chromosome(aln), labels=labels)
    checkIdentical(updt, chromosome(renew(aln, chromosome=updt)))

    obs <- renew(aln, chromosome=updt, position=1L+position(aln))
    checkIdentical(updt, chromosome(obs))
    checkIdentical(1L+position(aln), position(obs))

    checkException(renew(aln, position=1L), silent=TRUE)
}
