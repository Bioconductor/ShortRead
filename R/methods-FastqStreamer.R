.FastqStreamer_g$methods(
    add = function(bin) {
        if (verbose) {
            status(update=TRUE)
            msg("FastqStreamer$add()")
        }
        .Call(.streamer_add, sampler, bin, c(skips[ith], adds[ith]))
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

        if (!recycle && ith == length(skips))
            return (ShortReadQ())
        ith <<- ith %% length(skips) + 1L

        status(update=TRUE)

        prevTot <- status()["total"]
        if (status()["current"] != adds[ith]) {
            ## use C scratch buffer
            if (verbose)
                msg("FastqStreamer$yield() reader")
            add(raw())
        }

        while (0L != (adds[ith] - status()["current"])) {
            ## fill C buffer
            if (verbose)
                msg("FastqStreamer$yield() reader")
            bin <- reader(con, readerBlockSize)
            if (!length(bin))
                break

            currTot <- status()["total"]
            skips[ith] <<- max(0L, skips[ith] - (currTot - prevTot))
            prevTot <- currTot
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

setMethod(FastqStreamer, c("ANY", "missing"),
    function(con, n, readerBlockSize=1e8, verbose=FALSE)
{
    callGeneric(con, n=1e6, readerBlockSize=readerBlockSize,
        verbose=verbose)
})
          
setMethod(FastqStreamer, c("ANY", "numeric"),
    function(con, n, readerBlockSize=1e8, verbose=FALSE)
{
    n <- as.integer(n)
    if (length(n) != 1L || !is.finite(n) || n < 0L)
        stop("'n' must be length 1, finite and >= 0")
    if (is.character(con))
        con <- file(con)
    open(con, "rb")
    streamer <- .Call(.sampler_new, n)
    .ShortReadFile(.FastqStreamer_g, con, reader=.binReader,
        readerBlockSize=as.integer(readerBlockSize),
        skips = 0L, adds = n, ith = 0L, recycle=TRUE,
        sampler=streamer, verbose=verbose)
})

setMethod(FastqStreamer, c("ANY", "IRanges"),
    function(con, n, readerBlockSize=1e8, verbose=FALSE)
{
    if (is.character(con))
        con <- file(con)
    open(con, "rb")

    skips <- start(n) - c(1L, end(n)[-length(n)] + 1L)
    if (any(skips < 0)) {
        close(con)
        msg <- "'n' must have all(start(n)[-1] > end(n)[-length(n)])"
        .throw(SRError("UserArgumentMismatch", msg))
    }
    adds <- width(n)
    if (any(adds == 0)) {
        close(con)
        msg <- "'n' must have non-zero width()"
        .throw(SRError("UserArgumentMismatch", msg))
    }

    streamer <- .Call(.sampler_new, max(adds))
    .ShortReadFile(.FastqStreamer_g, con, reader=.binReader,
        readerBlockSize=as.integer(readerBlockSize),
        skips = skips, adds = adds, ith = 0L, recycle = FALSE, 
        sampler=streamer, verbose=verbose)
})

setMethod("FastqStreamerList", "ANY", 
          function(..., n, readerBlockSize=1e8, verbose=FALSE)
{
    FastqFileList(..., class="FastqStreamer")
})

setMethod("FastqStreamerList", "character",
          function(..., n, readerBlockSize=1e8, verbose=FALSE)
{
    listData <-
        lapply(..1, FastqStreamer, n=n, readerBlockSize=readerBlockSize,
               verbose=verbose)
    new("FastqStreamerList", listData=listData)
})
