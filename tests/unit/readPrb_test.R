sp <- SolexaPath(system.file("extdata", package="ShortRead"))

test_readPrb_input <- function()
{
    check <- function(obj) {
        checkTrue(validObject(obj))
        checkIdentical(256L, length(obj))
        checkIdentical(36L, unique(width(obj)))
    }
    icheck <- function(obj) {
        checkTrue(validObject(obj))
        checkIdentical(c(256L, 36L), dim(obj))
    }
    acheck <- function(obj, width) {
        checkIdentical("array", class(obj))
        checkIdentical("integer", typeof(obj))
        checkIdentical(c(256L, 4L, width), dim(obj))
    }
    check(readPrb(sp, ".*prb.txt", as="SolexaEncoding"))
    check(readPrb(sp, ".*prb.txt", as="FastqEncoding"))
    icheck(readPrb(sp, ".*prb.txt", as="IntegerEncoding"))
    acheck(readPrb(sp, ".*prb.txt", as="array"), 36L)
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
