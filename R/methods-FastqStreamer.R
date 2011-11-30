.FastqStreamer_g$methods(
    add = function(bin) {
        if (verbose) {
            status(update=TRUE)
            msg("FastqStreamer$add()")
        }
        .Call(.streamer_add, sampler, bin)
        status(update=TRUE)
    },
    status = function(update=FALSE) {
        "report status of FastqSampler"
        if (update || !length(.status))
            .status <<- .Call(.streamer_status, sampler)
        .status
    },
    yield = function(...) {
        "read at most n records in a connection"
        if (verbose) msg("FastqStreamer$yield()")

        status(update=TRUE)
        if (status()["current"] != status()["n"]) {
            ## use C scratch buffer
            if (verbose) {
                status(update=TRUE)
                msg("FastqStreamer$yield() reader")
            }
            add(raw())
        }

        while (status()["current"] != status()["n"]) {
            ## fill C buffer
            if (verbose) {
                status(update=TRUE)
                msg("FastqStreamer$yield() reader")
            }
            bin <- reader(con, readerBlockSize)
            if (!length(bin))
                break
            add(bin)
        }

        if (verbose)
            msg("FastqStreamer$yield() XStringSet")
        elts <- .Call(.streamer_as_XStringSet, sampler)
        if (verbose)
            msg("FastqStreamer$yield() ShortReadQ")
        ShortReadQ(elts[["sread"]], elts[["quality"]], elts[["id"]],
                   ...)
    })

FastqStreamer <-
    function(con, n=1e6, readerBlockSize=1e8, verbose=FALSE)
{
    if (length(n) != 1 || !is.finite(n) || n < 0)
        stop("'n' must be finite and >= 0")
    if (is.character(con))
        con <- file(con)
    open(con, "rb")
    streamer <- .Call(.sampler_new, as.integer(n))
    .FastqStreamer_g$new(con=con, reader=.binReader,
                         readerBlockSize=as.integer(readerBlockSize),
                         sampler=streamer, verbose=verbose)
}
