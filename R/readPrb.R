.readPrb <- function(file, ..., verbose=FALSE)
{
    if (verbose)
        cat(".prbScore", file, "\n")
    tbl <- read.table(file, colClasses="integer")
    cycles <- ncol(tbl) / 4
    if (!all.equal(cycles, as.integer(cycles)))
        .throw(SRError("Input/Output",
                       sprintf("%s\n  file: %s\n  columns: %d",
                               "'prb' column number not divisible by 4",
                               file, ncol(tbl))))
    score <- sapply(seq_len(cycles)-1, function(cyc, tbl) {
        do.call(pmax, tbl[, cyc*4+1:4])
    }, tbl)
    apply(matrix(as.raw(score+64), ncol=cycles),
          1, rawToChar)
}

.readPrb_character <- function(dirPath, pattern, ...)
{
    fls <- list.files(dirPath, pattern, full=TRUE)
    SFastqQuality(unlist(srapply(fls, .readPrb)))
}

setMethod("readPrb", "character", .readPrb_character)
