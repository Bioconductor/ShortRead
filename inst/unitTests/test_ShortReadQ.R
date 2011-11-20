sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")

checkShortReadQ <- function(obj, len, wd) {
    checkStringSet <- function(obj, type, len, wd) {
        checkTrue(is(obj, type))
        checkEquals(len, length(obj))
        checkEquals(wd, unique(width(obj)))
    }
    checkStringSet(obj, "ShortReadQ", len, wd[[1]])
    checkStringSet(sread(obj), "DNAStringSet", len, wd[[2]])
    checkStringSet(id(obj), "BStringSet", len, wd[[3]]) # ids w/ diff lengths
    checkStringSet(quality(obj), "QualityScore", len, wd[[4]])
}

.equals <- function(x, y)
{
    checkIdentical(as.character(sread(x)), as.character(sread(y)))
    checkIdentical(as.character(quality(quality(x))),
                   as.character(quality(quality(y))))
    checkIdentical(as.character(id(x)), as.character(id(y)))
}

test_ShortReadQ_constructors <- function() {
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    sr <- obj <- readFastq(sp)
    checkTrue(validObject(obj))
    checkShortReadQ(obj, 256, list(36, 36, 24:22, 36))

    obj <- ShortReadQ()
    checkTrue(class(obj) == "ShortReadQ")
    checkTrue(validObject(obj))

    obj <- ShortReadQ(sread(sr), quality(sr))
    checkTrue(class(obj) == "ShortReadQ")
    checkTrue(validObject(obj))
    .equals(new("ShortReadQ", sread=sread(sr),
                id=BStringSet(rep("", length(sr))),
                quality=quality(sr)), obj)

    obj <- ShortReadQ(sread(sr), quality(sr), id(sr))
    checkTrue(class(obj) == "ShortReadQ")
    checkTrue(validObject(obj))
    .equals(sr, obj)
}

test_FastqSampler_downsample <- function()
{
    sampler <- FastqSampler(file())
    downsample <- sampler$.downsample

    n <- length(letters)
    lets <- as.list(letters)
    LETS <- as.list(LETTERS)

    sampler$n <- n; sampler$tot_n <- 0L
    checkIdentical(lets, downsample(list(), lets))

    sampler$n <- sampler$tot_n <- n
    checkIdentical(lets, downsample(list(), lets))
    checkIdentical(lets, downsample(lets, list()))

    sampler$n <- 10L; sampler$tot_n <- n
    ans <- downsample(list(), lets)
    checkIdentical(sampler$n, length(ans))

    sampler$n <- sampler$tot_n <- 2L * n
    checkIdentical(c(lets, LETS), downsample(lets, LETS))

    sampler$n <- 2L * n + 1L
    checkIdentical(c(lets, LETS), downsample(lets, LETS))

    set.seed(123L)
    sampler$tot_n <- 2L * n; sampler$n <- n
    ans <- unlist(downsample(lets, LETS))
    checkIdentical(n, length(ans))
    checkIdentical(14L, sum(ans %in% letters))

    sampler$tot_n <- 200L * n
    ans <- unlist(downsample(lets, LETS))
    checkIdentical(n, length(ans))
    checkIdentical(n, sum(ans %in% lets))

    ans <- unlist(downsample(LETS, lets))
    checkIdentical(n, length(ans))
    checkIdentical(n, sum(ans %in% LETS))
}

test_FastqSampler <- function()
{
    sr <- readFastq(fl)
    ## here to re-use equality checker
    obj <- yield(FastqSampler(fl))
    .equals(sr, obj)

    yld <- yield(FastqSampler(fl, readerBlockSize=1000))
    checkTrue(validObject(yld))

    ## regression
    yld <- yield(FastqSampler(fl, readerBlockSize=256))
    checkIdentical(256L, length(yld))

}

test_FastqStreamer <- function()
{
    sr <- readFastq(fl)

    f <- FastqStreamer(fl, n=50)
    i <- 0L; len <- 0L
    while (length(y <- yield(f))) {
        len <- len + length(y)
        i <- i + 1L
    }
    checkIdentical(6L, i)
    checkIdentical(256L, len)

    ## values equal?
    f <- FastqStreamer(fl, n=50)
    .equals(sr[1:50], yield(f))
    .equals(sr[50+1:50], yield(f))

    ## whole file
    f <- FastqStreamer(fl, n=500)
    i <- 0L; len <- 0L
    while (length(y <- yield(f))) {
        .equals(sr, y)
        len <- len + length(y)
        i <- i + 1L
    }
    checkIdentical(1L, i)
    checkIdentical(256L, len)

    ## small reader block size
    f <- FastqStreamer(fl, n=50, readerBlockSize=100)
    i <- 0L; len <- 0L
    while (length(y <- yield(f))) {
        len <- len + length(y)
        i <- i + 1L
    }
    checkIdentical(6L, i)
    checkIdentical(256L, len)
}

test_ShortReadQ_coerce_QualityScaledDNAStringSet <- function()
{
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    obj <- readFastq(sp, qualityType="SFastqQuality")

    res <- as(obj, "QualityScaledDNAStringSet")
    checkIdentical(as.character(sread(obj)),
                   as.character(as(res, "DNAStringSet")))
    checkIdentical(as.character(quality(quality(obj))),
                   as.character(quality(res)))
    checkTrue(is(quality(res), "SolexaQuality"))

    obj <- initialize(obj, quality=FastqQuality(quality(quality(obj))))
    res <- as(obj, "QualityScaledDNAStringSet")
    checkIdentical(as.character(sread(obj)),
                   as.character(as(res, "DNAStringSet")))
    checkIdentical(as.character(quality(quality(obj))),
                   as.character(quality(res)))
    checkTrue(is(quality(res), "PhredQuality"))

    q <- MatrixQuality(as(quality(obj), "matrix"))
    obj <- initialize(obj, quality=q)
    checkException(as(obj, "QualityScaledDNAStringSet"), silent=TRUE)
}

test_ShortReadQ_subset <- function() {
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    obj <- readFastq(sp)

    obj1 <- obj[c(3, 5:7, 9)]
    checkShortReadQ(obj1, 5, list(36, 36, 23, 36))

    checkException(obj[,1], silent=TRUE)
    checkException(obj[1,1], silent=TRUE)
    checkIdentical(2L, length(obj[1:2,]))
    checkIdentical(2L, length(obj[1:2,drop=TRUE]))
    checkIdentical(2L, length(obj[1:2,,drop=TRUE]))
}

test_ShortReadQ_narrow <- function() {
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    sr <- readFastq(sp)

    obj <- narrow(sr, start=1, end=10)
    checkTrue(class(obj) == "ShortReadQ")
    checkTrue(length(obj) == length(sr))
    checkTrue(unique(width(obj)) == 10)
    checkIdentical(as.character(sread(obj)),
                   substr(as.character(sread(sr)), 1, 10))
    checkIdentical(as.character(quality(quality(obj))),
                   substr(as.character(quality(quality(sr))), 1, 10))
    checkIdentical(as.character(id(obj)), as.character(id(sr)))

    checkIdentical(narrow(sr, start=start(sread(sr))), sr)
}

test_ShortReadQ_compact <- function() {
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    sr <- readFastq(sp)[1:10]
    res <- compact(sr)
    checkIdentical(as.character(sread(sr)), as.character(sread(res)))
    checkIdentical(as.character(quality(quality(sr))),
                   as.character(quality(quality(res))))
}

test_ShortReadQ_clean <- function() {
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    obj <- readFastq(sp)

    cln <- clean(obj)
    checkIdentical(class(obj), class(cln))
    ## FIXME: need a stronger test
    checkEquals(length(obj), length(clean(obj)))
}

test_ShortReadQ_srsort <- function() {
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    obj <- readFastq(sp)
    srt <- srsort(obj)
    checkIdentical(class(obj), class(srt))
    checkIdentical(length(obj), length(srt))
    checkIdentical(srsort(sread(obj)), sread(srt))
    checkIdentical(quality(obj)[srorder(obj)], quality(srt))
}
