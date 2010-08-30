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
