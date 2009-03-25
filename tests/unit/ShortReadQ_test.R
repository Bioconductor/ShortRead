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

test_ShortReadQ_subset <- function() {
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    obj <- readFastq(sp)

    obj1 <- obj[c(3, 5:7, 9)]
    checkShortReadQ(obj1, 5, list(36, 36, 23, 36))

    checkException(obj[,1], silent=TRUE)
    checkException(obj[1,1], silent=TRUE)
    checkException(obj[1,], silent=TRUE)
    checkException(obj[1,], silent=TRUE)
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
