.readPrb <- function(file, ..., verbose=FALSE)
{
    if (verbose)
        cat(".prbScore", file, "\n")
    ln <- readLines(file, 1)
    cycles <- length(gregexpr("\t", ln, fixed=TRUE)[[1]]) + 1L
    tryCatch({
        .Call(.read_prb_as_character, file, cycles)
    }, error=function(err) {
        .throw(SRError("Input/Output",
                       sprintf("parsing 'prb'\n  file: %s\n  error: %s",
                               file,
                               conditionMessage(err))))
    })
}

.readPrb_character <- function(dirPath, pattern, ...)
{
    fls <- list.files(dirPath, pattern, full=TRUE)
    SFastqQuality(unlist(srapply(fls, .readPrb)))
}

setMethod("readPrb", "character", .readPrb_character)
