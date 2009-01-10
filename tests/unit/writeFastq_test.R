sp <- SolexaPath(system.file("extdata", package="ShortRead"))

test_writeFastq_roundtrip <- function()
{
    ## potential coercion from '.' to 'N'
    rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
    file <- tempfile()
    writeFastq(rfq, file)
    fq <- readFastq(dirname(file), basename(file))

    checkIdentical(as.character(id(fq)), as.character(id(rfq)))
    checkIdentical(as.character(sread(fq)), as.character(sread(rfq)))
    checkIdentical(as.character(quality(quality(fq))),
                   as.character(quality(quality(rfq))))
}
