## Function: append (package base)
## x="AlignedDataFrame", values="AlignedDataFrame", after="missing"
## x="AlignedRead", values="AlignedRead", after="missing"
## x="QualityScore", values="QualityScore", after="missing"
## x="ShortRead", values="ShortRead", after="missing"
## x="ShortReadQ", values="ShortReadQ", after="missing"

sp <- SolexaPath(system.file("extdata", package="ShortRead"))

.equal <- function(x, y)
{
    checkIdentical(class(x), class(y))
    checkIdentical(length(x), length(y))
    checkIdentical(as.character(id(x)), as.character(id(y)))
    checkIdentical(as.character(sread(x)), as.character(sread(y)))
    if (is(x, "ShortReadQ"))
        checkIdentical(as.character(quality(quality(x))),
                       as.character(quality(quality(y))))
    if (is(x, "AlignedRead")) {
        checkIdentical(strand(x), strand(y))
        checkIdentical(chromosome(x), chromosome(y))
        checkIdentical(position(x), position(y))
        checkIdentical(dim(alignData(x)), dim(alignData(y)))
        adx <- alignData(x); ady <- alignData(y)
        checkIdentical(varMetadata(adx), varMetadata(ady))
        pdx <- pData(adx); pdy <- pData(ady)
        row.names(pdx) <- row.names(pdy) <- NULL
        checkIdentical(pdx, pdy)
    }
}

test_append <- function() 
{
    aln <- readAligned(sp, "s_2_export.txt")
    aaln <- append(aln, aln)
    .equal(aln, aaln[seq_len(length(aln))])
    .equal(aln, aaln[length(aln) + seq_len(length(aln))])
}

test_append_exception <- function()
{
    checkException(append(sp, sp), silent=TRUE)
    aln <- readAligned(sp, "s_2_export.txt")
    checkException(append(quality(aln), aln), silent=TRUE)
}
