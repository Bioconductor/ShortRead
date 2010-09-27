## see test_ShortReadQ

test_Sampler_rec_counter <- function()
{
    CFUN <- ShortRead:::.sampler_rec_counter

    sep <- charToRaw("@")
    checkIdentical(3L, .Call(CFUN, charToRaw("@one@two@three"), sep))
    checkIdentical(3L, .Call(CFUN, charToRaw("@one@two@"), sep))
    checkIdentical(3L, .Call(CFUN, charToRaw("@@@"), sep))
    ## first record is always tagged, regardless of delimiter
    checkIdentical(3L, .Call(CFUN, charToRaw("one@two@three"), sep))
    sep <- charToRaw("\n@")
    checkIdentical(3L, .Call(CFUN, charToRaw("@one\n@two\n@three"), sep))
    checkIdentical(2L, .Call(CFUN, charToRaw("@one\ntwo\n@three"), sep))
    checkIdentical(2L, .Call(CFUN, charToRaw("one\ntwo\n@three"), sep))
}

test_Sampler_rec_parser <- function()
{
    CFUN <- ShortRead:::.sampler_rec_parser

    char <- "@one@two@three"
    buf <- charToRaw(char)
    sep <- charToRaw("@")
    exp <- lapply(list("@one", "@two", "@three"), charToRaw)
    checkIdentical(exp, .Call(CFUN, buf, sep, 1:3))
    checkIdentical(exp[1], .Call(CFUN, buf, sep, 1L))
    checkIdentical(exp[c(1,3)], .Call(CFUN, buf, sep, c(1L, 3L)))
    checkIdentical(exp[3], .Call(CFUN, buf, sep, 3L))
}

## test_Sampler_as_fastsq <- function()
## {
##     CFUN <- ShortRead:::.sampler_as_fastq
## }
