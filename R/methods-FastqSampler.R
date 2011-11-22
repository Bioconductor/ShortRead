.FastqSampler_g$methods(
    initialize = function(..., n) {
        callSuper(...)
        sampler <<- .Call(.sampler_new, as.integer(n))
        .self
    },
    reset = function() {
        "reopen the connection"
        if (verbose) msg("FastqSampler$reset()")
        if (isOpen(con)) {
            if (verbose) msg("FastqSamper$reset() re-open")
            s <- summary(con)
            class <- s$class
            desc <- s$description
            close(con)
            con <<- do.call(s$class, list(desc, "rb"))
      } else {
          open(con, "rb")
      }
      .self
    },
    yield=function() {
        "read and sample all records in a connection"
        if (verbose) msg("FastqSampler$yield()")
        reset()
        while (length(bin <- reader(con, readerBlockSize))) {
            if (verbose) {
                status(update=TRUE)
                msg("FastqSampler$yield() reader")
            }
            .Call(.sampler_add, sampler, bin)
        }
        status(update=TRUE)
        if (verbose)
            msg("FastqSampler$yield() as XStringSet")
        elts <- .Call(.sampler_as_XStringSet, sampler)
        ShortReadQ(elts[["sread"]], elts[["quality"]], elts[["id"]])
    },
    status=function(update=FALSE) {
        "report status of FastqSampler"
        if (update || !length(status_txt))
            status_txt <<- .Call(.sampler_summary, sampler)
        status_txt
    },
    show = function() {
        cat("class:", class(.self), "\n")
        cat("file:", basename(summary(.self$con)$description), "\n")
        s <- .self$status()
        cat("status:", paste(names(s), s, sep="=", collapse=" "), "\n")
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
                        verbose=verbose)
}

setMethod(yield, "FastqSampler", function(x, ...) x$yield())
