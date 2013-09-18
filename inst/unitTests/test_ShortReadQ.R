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

test_FastqSampler <- function()
{
    sr <- readFastq(fl)
    ## here to re-use equality checker
    
    obj <- yield(f <- FastqSampler(fl))
    close(f)
    .equals(sr, obj)

    yld <- yield(f <- FastqSampler(fl, readerBlockSize=1000))
    close(f)
    checkTrue(validObject(yld))

    ## regression
    yld <- yield(f <- FastqSampler(fl, readerBlockSize=256))
    close(f)
    checkIdentical(256L, length(yld))

}

test_FastqSampler_rand <- function()
{
    ## two samples with the same random number seed are identical
    samp <- FastqSampler(fl, 50)
    set.seed(123L); obs <- yield(samp)
    set.seed(123L); exp <- yield(samp)
    close(samp)
    .equals(obs, exp)

    ## different samples
    set.seed(123L)
    samp <- FastqSampler(fl, 50)
    obs <- length(Reduce(intersect, replicate(2, id(yield(samp)))))
    checkIdentical(7L, obs)
    obs <- length(Reduce(intersect, replicate(3, id(yield(samp)))))
    checkIdentical(0L, obs)
    close(samp)
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
    close(f)
    checkIdentical(6L, i)
    checkIdentical(256L, len)

    ## values equal?
    f <- FastqStreamer(fl, n=50)
    .equals(sr[1:50], yield(f))
    .equals(sr[50+1:50], yield(f))
    close(f)

    ## whole file
    f <- FastqStreamer(fl, n=500)
    i <- 0L; len <- 0L
    while (length(y <- yield(f))) {
        .equals(sr, y)
        len <- len + length(y)
        i <- i + 1L
    }
    close(f)
    checkIdentical(1L, i)
    checkIdentical(256L, len)

    ## small reader block size
    f <- FastqStreamer(fl, n=50, readerBlockSize=100)
    i <- 0L; len <- 0L
    while (length(y <- yield(f))) {
        len <- len + length(y)
        i <- i + 1L
    }
    close(f)
    checkIdentical(6L, i)
    checkIdentical(256L, len)
}

test_FastqStreamer_roundtrip <- function()
{
    out <- tempfile()
    writeFastq(v1 <- readFastq(fl), out)
    s <- FastqStreamer(out)
    .equals(v1, yield(s))
}

test_FastqStreamer_IRanges <- function()
{
    sr <- readFastq(fl)

    ## basics
    rng <- IRanges(c(50, 100, 200), width=c(5, 4, 3))
    f <- FastqStreamer(fl, rng)
    .equals(sr[50:54], yield(f))
    .equals(sr[100:103], yield(f))
    .equals(sr[200:202], yield(f))
    .equals(ShortReadQ(), yield(f))
    close(f)

    ## successive
    rng <- IRanges(c(50, 60), width=10)
    f <- FastqStreamer(fl, rng)
    .equals(sr[50:59], yield(f))
    .equals(sr[60:69], yield(f))
    .equals(ShortReadQ(), yield(f))
    close(f)

    ## off-the-end
    rng <- IRanges(250, width=100)
    f <- FastqStreamer(fl, rng)
    .equals(sr[250:256], yield(f))
    .equals(ShortReadQ(), yield(f))
    close(f)

    ## too-short buffer to skip all reads in one binary input
    rng <- IRanges(250, width=5)
    f <- FastqStreamer(fl, rng, readerBlockSize=10000)
    .equals(sr[250:254], yield(f))
    .equals(ShortReadQ(), yield(f))
    close(f)

    rng <- IRanges(241, width=5)
    f <- FastqStreamer(fl, rng, readerBlockSize=10000)
    .equals(sr[241:245], yield(f))
    .equals(ShortReadQ(), yield(f))
    close(f)

    ## exceptions
    rng <- IRanges(50, 49)              # non-zero
    checkException(FastqStreamer(fl, rng), silent=TRUE)
    rng <- IRanges(c(50, 59), c(60, 70)) # strictly increasing
    checkException(FastqStreamer(fl, rng), silent=TRUE)
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

test_ShortReadQ_coerce_matrix <- function() 
{
    ## 0-length
    fq <- FastqQuality()
    exp <- matrix(NA_integer_, 0, 0)
    checkIdentical(exp, as(fq, "matrix"))

    ## ragged matrix
    fq <- FastqQuality(BStringSet(c("]]X", "]]]X")))
    exp <- matrix(c(rep(60L, 4), 55L, 60L, NA_integer_, 55L), 2)
    checkIdentical(exp, as(fq, "matrix"))
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
