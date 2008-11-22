.readPrb <- function(file, ..., asSolexa, verbose)
{
    if (verbose)
        cat(".readPrb", file, "\n")
    tryCatch({
        .Call(.read_prb_as_character, file, asSolexa)
    }, error=function(err) {
        .throw(SRError("Input/Output",
                       sprintf("parsing 'prb'\n  file: %s\n  error: %s",
                               file,
                               conditionMessage(err))))
    })
}

.readPrb_quality <-
    function(dirPath, pattern, qclass, ..., asSolexa, verbose)
{
    fls <- .file_names(dirPath, pattern)
    qclass(unlist(srapply(fls, .readPrb, ...,
                          asSolexa=asSolexa, verbose=verbose)))
}

.readPrb_IntegerEncoding <-
    function(dirPath, pattern, ..., verbose)
{
    res <- .readPrb_quality(dirPath, pattern, SFastqQuality, ...,
                            asSolexa=TRUE, verbose=verbose)
    if (length(unique(width(res)))!=1)
        .throw(SRError("Input/Output", "reads have different widths") )
    as(res, "matrix")
}

.readPrb_array <-
    function(dirPath, pattern, ..., verbose=FALSE)
{
    nrec <- countLines(dirPath, pattern)
    crec <- c(0, cumsum(nrec))
    fls <- .file_names(dirPath, pattern)
    gz <- gzfile(fls[[1]], "rb")
    tryCatch({
        ln <- readLines(gz, 1)
    }, finally=close(gz))
    cycles <- length(gregexpr("\t", ln, fixed=TRUE)[[1]]) + 1L    
    a <- array(integer(), c(sum(nrec), 4L, cycles),
               dimnames=list(NULL, c("A", "C", "G", "T"), NULL))
    what <- rep(list(integer()), 4L * cycles)
    for (i in seq_along(fls))
        tryCatch({
            gz <- gzfile(fls[[i]], "rb")
            data <- unlist(scan(gz, what, sum(nrec), ..., quiet=!verbose))
            a[(crec[i]+1):crec[i+1],,] <-
                array(data, c(nrec[[i]], 4L, cycles))
        }, error=function(err) {
            .throw(SRError("Input/Output",
                           sprintf("parsing 'prb'\n  file: %s\n  error: %s",
                                   fls[[i]],
                                   conditionMessage(err))))
        }, finally=close(gz))
    a
}

.readPrb_character <-
    function(dirPath, pattern=character(0),
             as=c(
               "SolexaEncoding", "FastqEncoding", "IntegerEncoding",
               "array"),
             ..., verbose=FALSE)
{
    if (missing(as)) {
        as <- "SolexaEncoding"
    } else if (!is.character(as) || length(as) != 1) {
        .arg_mismatch_type_err("as", "character(1)")
    } else {
        vals <- eval(formals(.readPrb_character)$as)
        if (!as %in% vals)
            .arg_mismatch_value_err("as", as, vals)
    }
    tryCatch({
        switch(as,
               SolexaEncoding=.readPrb_quality(
                 dirPath, pattern, SFastqQuality, ..., asSolexa=TRUE,
                 verbose=verbose),
               FastqEncoding=.readPrb_quality(
                 dirPath, pattern, FastqQuality, ..., asSolexa=FALSE,
                 verbose=verbose),
               IntegerEncoding=.readPrb_IntegerEncoding(
                 dirPath, pattern, ..., verbose=verbose),
               array=.readPrb_array(
                 dirPath, pattern, ..., verbose=verbose))
    }, error=function(err) {
        if (is(err, "SRError")) stop(err)
        else {
            pat <- paste(pattern, collapse=" ")
            txt <- paste("'%s' failed to parse files",
                         "dirPath: '%s'",
                         "pattern: '%s'",
                         "as: '%s'",
                         "error: %s", sep="\n  ")
            msg <- sprintf(txt, "readPrb",
                           paste(dirPath, collapse="'\n    '"),
                           pat, as, conditionMessage(err))
            .throw(SRError("Input/Output", msg))
        }
    })
}

setMethod("readPrb", "character", .readPrb_character)
