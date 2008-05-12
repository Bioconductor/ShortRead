test_AlignedRead_readAligned_SolexaExport <- function() {
    fl <- system.file("extdata", "Data", "C1-36Firecrest", "Bustard",
                      "GERALD", package="ShortRead")
    obj <- readAligned(fl, pattern="s_2_export.txt", type="SolexaExport")
    checkTrue(validObject(obj))
    checkTrue(is(quality(obj), "SFastqQuality"))
    checkTrue(is(alignQuality(obj), "NumericQuality"))
    checkIdentical(varLabels(alignData(obj)),
                   c("run", "lane", "tile", "x", "y", "filtering"))
}

test_AlignedRead_readAligned_MAQMapview <- function() {
    fl <- system.file("extdata", "maq", package="ShortRead")
    obj <- readAligned(fl, pattern=".*aln.*", type="MAQMapview")
    checkTrue(validObject(obj))
    checkTrue(is(quality(obj), "FastqQuality"))
    checkTrue(is(alignQuality(obj), "NumericQuality"))
    checkIdentical(varLabels(alignData(obj)),
                   c("nMismatchBestHit", "mismatchQuality",
                   "nExactMatch24", "nOneMismatch24"))
}
