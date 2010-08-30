test_srapply_USE.NAMES <- function()
{
    f <- function(x, ...) x
    checkIdentical(NULL, names(srapply(letters, f, USE.NAMES=FALSE)))
    checkIdentical(letters,
                   names(srapply(letters, f, USE.NAMES=TRUE)))
}

test_srapply_nested_error <- function()
{
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    fl <- analysisPath(sp)
    suppressWarnings(tryCatch({
        qa(fl, "s_2_export.txt", "fastq")
    }, SRError=function(err) {
        txt <- "ValueUnavailable\n  0 elements returned; expected >=1"
        checkIdentical(txt, err$message)
    }))
}
