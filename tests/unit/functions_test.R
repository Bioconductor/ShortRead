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

## sort / order

test_alphabetOrder <- function() {
    ## setup
    oldc <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", oldc))
    Sys.setlocale("LC_COLLATE", "C")
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")

    checkEquals(srorder(sread(rfq)),
                order(as.character(sread(rfq))))
    checkEquals(srorder(quality(rfq)),
                order(as.character(quality(quality(rfq)))))

    checkEquals(srduplicated(sread(rfq)),
                duplicated(as.character(sread(rfq))))
    checkEquals(srduplicated(quality(rfq)),
                duplicated(as.character(quality(quality(rfq)))))
}
