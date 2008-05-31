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

readAligned_maq_consistent <- function() {
    ## FIXME: find adequate data to store in ShortRead pkg
    if (!file.exists("/home/jdavison/sharedrsrc/proj/ycao/data/binary_maps/s_5.map"))
        return(TRUE)

    x <- readAligned("/home/jdavison/sharedrsrc/proj/ycao/data/binary_maps",
                     "s_5.map", "MAQMap")
    y <- readAligned("/home/jdavison/sharedrsrc/proj/ycao/data/text_maps",
                     "s_5.txt", "MAQMapview")

    checkIdentical(length(x), length(y))
    checkIdentical(width(x), width(y))
    checkIdentical(as.character(chromosome(x)), as.character(chromosome(y)))

    ## FIXME: we'd really like chromosome to have identical levels,
    ## but info on levels with no mapped reads is not available in the
    ## text version
    idx <- match(levels(chromosome(y)), levels(chromosome(x)))
    checkTrue(all(!is.na(idx)) && all(diff(idx) > 0))
    checkIdentical(position(x), position(y))
    checkIdentical(strand(x), strand(y))
    checkIdentical(alignQuality(x), alignQuality(y))
    checkIdentical(alignData(x), alignData(y))

    .checkXString <- function(x, y) {
        checkIdentical(as.character(x), as.character(y))
    }
    .checkXString(sread(x), sread(y))
    .checkXString(sread(x), sread(y))
    .checkXString(quality(quality(x)), quality(quality(y)))
    .checkXString(id(x), id(y))
}
