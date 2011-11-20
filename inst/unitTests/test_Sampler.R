## see test_ShortReadQ

test_Sampler_rec_counter <- function()
{
    CFUN <- ShortRead:::.sampler_rec_counter

    str <- "@one\nAA\n+\nII\n@two\nAA\n+\nII\n@three\nAA\n+\nII"
    checkIdentical(3L, .Call(CFUN, charToRaw(str)))
    str <- "@one\nAA\n+\nII\n@two\nAA\n+\nII\n@three\nAA\n+\nII\n"
    checkIdentical(3L, .Call(CFUN, charToRaw(str)))
    ## counts partial records, too
    str <- "@one\nAA\n+\nII\n@two\nAA\n+\nII\n@three"
    checkIdentical(3L, .Call(CFUN, charToRaw(str)))
}

test_Sampler_rec_parser <- function()
{
    CFUN <- ShortRead:::.sampler_rec_parser

    str <- "@one\nAA\n+\nII\n@two\nAA\n+\nII\n@three\nAA\n+\nII"
    buf <- charToRaw(str)
    obs <- paste(sapply(.Call(CFUN, buf, 3L), rawToChar), collapse="")
    checkIdentical(str, obs)
}

## test_Sampler_as_fastsq <- function()
## {
##     CFUN <- ShortRead:::.sampler_as_fastq
## }
