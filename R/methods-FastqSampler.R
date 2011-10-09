.FastqSampler_g$methods(
    .add = function(bin, flush=FALSE) {
        ".add (incomplete) 'bin'ary stream, possibly flush'ing buffer"
        if (verbose) msg("FastqSampler$.add")
        res <- recParser(buf, bin, n, tot_n)
        samp <- res[["parsed_bin"]]
        if (flush) {
            buf <<- raw()
            if (tot_n > n && runif(1L) > 1 / tot_n) # sample buf?
                samp <- samp[-length(samp)]
        } else {
            len <- length(samp)
            buf <<- samp[[len]]
            res[["rec_n"]] <- res[["rec_n"]] - 1L
            samp <- samp[-len]
        }
        if (length(samp)) {
            if (length(records) + length(samp) <= n) {
                records <<- c(records, samp)
            } else if (length(records) < n) {
                len <- length(records) + length(samp) - n
                drop <- base::sample(length(records), len)
                records <<- c(records[-drop], samp)
            } else {
                len <- length(samp)
                records[base::sample(n, len)] <<- samp
            }
        }
        saved_n <<- length(records)
        tot_n <<- tot_n + res[["rec_n"]]
        .self
    },
    .yield = function() {
        "read and sample all records in a connection"
        if (verbose) msg("FastqSampler$.yield")
        reset()
        while (0 != length(bin <- reader(con, readerBlockSize)))
            .add(bin)
        flush()
        if (verbose) msg("FastqSampler$sample finished")
        .self
    })

FastqSampler <-
    function(con, n = 1e6, readerBlockSize=1e8, verbose=FALSE)
{
    if (length(n) != 1 || !is.finite(n) || n < 0)
        stop("'n' must be finite and >= 0")
    if (is.character(con))
        con <- file(con)
    .FastqSampler_g$new(con=con, n=as.integer(n), reader=.binReader,
                        readerBlockSize=as.integer(readerBlockSize),
                        recParser=.fixedBinRecSampler,
                        verbose=verbose)
}

setMethod(yield, "FastqSampler",
    function(x, ...)
{
    x$.yield()
    elts <- .Call(.sampler_as_fastq, x$get())
    ShortReadQ(elts[["sread"]], elts[["quality"]], elts[["id"]], ...)
})
