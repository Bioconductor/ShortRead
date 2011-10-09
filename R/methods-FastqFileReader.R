.binReader <-
    function(con, n)
    ## read 'n' bytes from 'con', returning raw()
{
    readBin(con, raw(), n)
}

.fixedBinRecSampler <-
    function(buf, bin, n, tot_n)
    ## sample c('buf', 'bin') records in proportion to their
    ## reprentation, returning an environment() of raw() parsed_bin
    ## records + new buf
{
    bin <- c(buf, bin)
    env <- new.env(parent=emptyenv())
    env[["rec_n"]] <- rec_n <- .Call(.sampler_rec_counter, bin)
    if (1L >= rec_n) {
        env[["parsed_bin"]] <- list(bin)
        return(env)
    }

    rec_n <- rec_n - 1L                 # buf as last element
    if (tot_n + rec_n <= n)             # all records
        samp <- seq_len(rec_n)
    else {                              # records in ppn to abundance
        trials <- min(n, rec_n)
        p <- rec_n / (tot_n + rec_n)
        samp <- sort(sample(rec_n, rbinom(1L, trials, p)))
    }
    env[["parsed_bin"]] <-
        .Call(.sampler_rec_parser, bin, c(samp, rec_n + 1L))
    env
}

