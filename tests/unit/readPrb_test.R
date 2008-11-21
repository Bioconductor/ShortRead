sp <- SolexaPath(system.file("extdata", package="ShortRead"))

test_readPrb_input <- function()
{
    check <- function(obj) {
        checkTrue(validObject(obj))
        checkIdentical(256L, length(obj))
        checkIdentical(36L, unique(width(obj)))
    }
    mcheck <- function(obj, width) {
        checkIdentical("matrix", class(obj))
        checkIdentical("integer", typeof(obj))
        checkIdentical(c(256L, width), dim(obj))
    }
    check(readPrb(sp, ".*prb.txt", as="SolexaEncoding"))
    check(readPrb(sp, ".*prb.txt", as="FastqEncoding"))
    mcheck(readPrb(sp, ".*prb.txt", as="IntegerEncoding"), 36L)
    mcheck(readPrb(sp, ".*prb.txt", as="matrix"), 144L)
}

test_readPrb_consistent <- function()
{
    exp <- readPrb(sp, ".*prb.txt", as="IntegerEncoding")
    checkIdentical(exp,
                   as(readPrb(sp, ".*prb.txt", as="SolexaEncoding"), "matrix"))
    checkIdentical(exp,
                   as(readPrb(sp, ".*prb.txt", as="FastqEncoding"), "matrix"))
}

test_readPrb_errors <- function()
{
}
