.FastqSampler_g$methods(
    .downsample = function(records, samp) {
        ##   |-- x = length(records) --|-- y = length(samp) --|
        ##   |---- n = length(result) ----|-------- w --------|
        x <- length(records); y <- length(samp)
        if (x + y <= n) {
            c(records, samp)
        } else {
            w <- x + y - n
            samp_n <- n - x + rbinom(1L, min(w, x), w / tot_n)
            if (samp_n) {
                records_n <- n - samp_n
                c(sample(records, records_n), sample(samp, samp_n))
            } else {
                records
            }
        }
    },
    .add = function(bin, flush=FALSE) {
        ".add (incomplete) 'bin'ary stream, possibly flush'ing buffer"
        if (verbose) msg("FastqSampler$.add")
        res <- recParser(buf, bin, n)
        samp <- res[["parsed_bin"]]
        tot_n <<- tot_n + res[["rec_n"]]
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
        records <<- .downsample(records, samp)
        saved_n <<- length(records)
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
                        recParser=.fixedBinRecParser,
                        verbose=verbose)
}

setMethod(yield, "FastqSampler",
    function(x, ...)
{
    recs <- x$.yield()$get()
    if (0L == length(recs))
        return(ShortReadQ())

    ## split many recs into smaller units to be append'ed (cheap)
    asShortReadQ <- function(idx, recs) {
        elts <- .Call(.sampler_as_fastq, recs[idx])
        ShortReadQ(elts[["sread"]], elts[["quality"]], elts[["id"]])
    }
    .binsize <- .Machine$integer.max
    bin <- floor(cumsum(as.numeric(sapply(recs, length))) / .binsize)
    idx <- split(seq_along(recs), bin)
    res <- asShortReadQ(idx[[1]], recs)
    for (i in idx[-1])
        res <- append(res, asShortReadQ(i, recs))
    res
})
