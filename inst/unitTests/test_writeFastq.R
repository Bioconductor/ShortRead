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

test_writeFastq_writeError <- function()
{
    object <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
    file <- tempfile()
    mode <- "w"
    max_width <- 10L
    .write_fastq <- ShortRead:::.write_fastq
    checkException(.Call(.write_fastq, id(object), sread(object),
                         quality(quality(object)), file, mode,
                         max_width),
                   silent=TRUE)
}

test_writeFastq_roundtrip0length <- function()
{
    dest <- tempfile()
    file.create(dest)
    exp <- readFastq(dest)
    writeFastq(exp, dest <- tempfile())
    checkIdentical(exp, readFastq(dest))
}
