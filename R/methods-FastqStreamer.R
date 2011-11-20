.FastqStreamer_g$methods(
    .add = function(bin, flush=FALSE) {
        ".add (incomplete) 'bin'ary stream, possibly flush'ing buffer"
        if (verbose) msg("FastqSampler$.add")
        res <- recParser(buf, bin, Inf)
        samp <- res[["parsed_bin"]]
        if (flush) {
            buf <<- raw()
        } else {
            len <- length(samp)
            buf <<- samp[[len]]
            res[["rec_n"]] <- res[["rec_n"]] - 1L
            samp <- samp[-len]
        }
        if (length(samp))
            records <<- c(records, samp)
        saved_n <<- length(records)
        tot_n <<- tot_n + res[["rec_n"]]
        .self
    },
    get = function() {
        "at most 'n' complete records"
        if (verbose) msg("FastqStreamer$get")

        while (n > length(records) &&
               0 != length(bin <- reader(con, readerBlockSize)))
            .add(bin)
        if (n > length(records))
            flush()                  # last record
        if (n < length(records)) {
            idx <- 1:n
            res <- records[idx]
            records <<- records[-idx]
        } else {
            res <- records
            records <<- list()
        }
        res
    })

FastqStreamer <-
    function(con, n=1e6, readerBlockSize=1e8, verbose=FALSE)
{
    if (length(n) != 1 || !is.finite(n) || n < 0)
        stop("'n' must be finite and >= 0")
    if (is.character(con))
        con <- file(con, "rb")
    .FastqStreamer_g$new(con=con, n=as.integer(n), tot_n=0L,
                         saved_n=0L, reader=.binReader,
                         readerBlockSize=as.integer(readerBlockSize),
                         recParser=.fixedBinRecParser,
                         verbose=verbose)
}

setMethod(yield, "FastqStreamer",
    function(x, ...)
{
    elts <- .Call(.sampler_as_fastq, x$get())
    if (length(elts[["sread"]]))
        ShortReadQ(elts[["sread"]], elts[["quality"]], elts[["id"]], ...)
    else
        ShortReadQ()
})
