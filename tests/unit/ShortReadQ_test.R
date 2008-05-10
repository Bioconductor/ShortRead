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

test_ShortReadQ_constructors <- function() {
    ## no direct construction

    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    obj <- readFastq(sp)
    checkTrue(validObject(obj))
    checkShortReadQ(obj, 256, list(36, 36, 24:22, 36))
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
