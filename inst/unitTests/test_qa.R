test_missingLaneName <- function()
{
    caught <- FALSE
    tryCatch(qa(AlignedRead()),
             error=function(err) {
                 caught <<- conditionMessage(err) == 
                     "UserArgumentMismatch\n  'lane' must be 'character(1)'"
             })
    checkTrue(caught)
}

test_no_replicate_reads <- function()
{
    df <- data.frame(nOccurrences=1, nReads=10, lane=1)
    x <- ShortRead:::.plotReadOccurrences(df)
    checkTrue(is(x, "trellis"))
}

test_qa_alphabetFrequency <- function()
{
    FUN <- ShortRead:::.qa_alphabetFrequency

    checkException(FUN(DNAStringSet()), silent=TRUE)

    exp <- alphabetFrequency(DNAStringSet(), collapse=TRUE,
                             baseOnly=TRUE)
    checkEquals(exp, FUN(DNAStringSet(), collapse=TRUE, baseOnly=TRUE))

    exp <- alphabetFrequency(DNAStringSet(), collapse=TRUE)
    checkEquals(exp,
                FUN(DNAStringSet(), collapse=TRUE))

    dna <- DNAStringSet(c("ACTG", "GTCANM"))
    checkEquals(alphabetFrequency(dna, collapse=TRUE),
                FUN(dna, collapse=TRUE))
    checkEquals(alphabetFrequency(dna, collapse=TRUE, baseOnly=TRUE),
                FUN(dna, collapse=TRUE, baseOnly=TRUE))
}
