.binReader <-
    function(con, n)
    ## read 'n' bytes from 'con', returning raw()
{
    readBin(con, raw(), n)
}

.fixedBinRecParser <-
    function(buf, bin, n)
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
    env[["parsed_bin"]] <- .Call(.sampler_rec_parser, bin, rec_n)
    env
}

