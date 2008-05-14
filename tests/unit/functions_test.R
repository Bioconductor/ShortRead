## readFastq

test_readFastq_errors <- function() {
    checkTrue(FALSE)
}

## alphabetByCycle

checkAlphabetByCycle <- function(obj) {
    abc <- alphabetByCycle(obj)
    validObject(abc)
    checkEquals(length(obj)*unique(width(obj)), sum(abc))
}

test_alphabetByCycle <- function() {
    sp <- SolexaPath(system.file('extdata', package="ShortRead"))
    sq <- readFastq(sp)

    checkAlphabetByCycle(sread(sq))
    checkAlphabetByCycle(quality(sq))

    obj <- alphabetByCycle(sq)
    validObject(obj)
    checkEquals(2, length(obj))
    checkEquals(alphabetByCycle(sread(sq)), obj[["sread"]])
    checkEquals(alphabetByCycle(quality(sq)), obj[["quality"]])

    checkException(alphabetByCycle(id(sq)), silent=TRUE) # unequal widths
    checkException(alphabetByCycle(sread(new("ShortReadQ"))),
                   silent=TRUE)         # zero-length
    checkException(alphabetByCycle(new("ShortReadQ")), silent=TRUE)
}

## countLines

test_countLines <- function() {
    sp <- SolexaPath(system.file('extdata', package="ShortRead"))
    nlines <- countLines(analysisPath(sp), "s_1_sequence.txt")
    exp <- 1024; names(exp) <- "s_1_sequence.txt"
    checkEquals(exp, nlines)
    checkException(countLines(tempdir()), silent=TRUE)
}
