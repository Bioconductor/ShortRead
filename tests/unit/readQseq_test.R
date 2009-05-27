sp <- SolexaPath(system.file("extdata", package="ShortRead"))

checkBstring <- function(obs, exp)
{
    checkEquals(as.character(obs), as.character(exp))
}

test_readQseq_ShortReadQ<- function()
{
    res <- readQseq(sp)
    checkEquals("ShortReadQ", as.vector(class(res)))
    checkEquals(256L, length(res))
    checkEquals(158223L, sum(alphabetScore(res)))
    alf <- alphabetFrequency(sread(res), collapse=TRUE, baseOnly=TRUE)
    checkEquals(structure(c(1697L, 1639L, 1481L, 1706L, 133L),
                          .Names = c("A", "C", "G", "T", "other")),
                alf)
}

test_readQseq_ShortReadQ_filtered <- function()
{
    res <- readQseq(sp, filtered=TRUE)
    checkEquals(187L, length(res))
}

test_readQseq_DataFrame <- function()
{
    res <- readQseq(sp)
    xdf <- readQseq(sp, as="DataFrame")
    checkEquals("DataFrame", as.vector(class(xdf)))
    checkEquals(c(256L, 11L), dim(xdf))
    checkBstring(sread(res), xdf[[9]])
    checkBstring(quality(quality(res)), xdf[[10]])
}

test_readQseq_DataFrame_filtered <- function()
{
    xdf0 <- readQseq(sp, as="DataFrame")
    xdf0 <- xdf0[xdf0[[11]]=="Y", -11]
    xdf <- readQseq(sp, as="DataFrame", filtered=TRUE)
    checkEquals(dim(xdf0), dim(xdf))
    for (i in 1:8)
        checkEquals(xdf0[[i]], xdf[[i]])
    for (i in 9:10)
        checkBstring(xdf0[[i]], xdf[[i]])
}
