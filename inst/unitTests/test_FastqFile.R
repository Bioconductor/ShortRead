sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")

test_FastqFile <- function()
{
    fq <- FastqFile(fl)
    checkTrue(validObject(fq))
    checkIdentical(path(fq), fl)
    checkTrue(!isOpen(fq))
    close(fq)
}

test_FastqFileList <- function()
{
    fql0 <- FastqFileList(c(fl, fl))
    checkTrue(validObject(fql0))
    checkIdentical(2L, length(fql0))
    fql1 <- FastqFileList(FastqFile(fl), FastqFile(fl))
    checkIdentical(sapply(fql0, path), sapply(fql1, path))
    checkIdentical(sapply(fql0, isOpen), sapply(fql1, isOpen))
    open(fql0)
    checkTrue(all(sapply(fql0, isOpen)))
    close(fql0)
    checkTrue(all(!sapply(fql0, isOpen)))
    checkTrue(all(!sapply(fql1, isOpen)))
}
