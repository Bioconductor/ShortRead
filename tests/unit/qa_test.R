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

test_emptyLane <- function()
{
    checkTrue(FALSE)
##     checkTrue(validObject(report(qa(AlignedRead(), "foo"))))
}
