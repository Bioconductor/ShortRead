.binReader <-
    function(con, n)
    ## read 'n' bytes from 'con', returning raw()
{
    readBin(con, raw(), n)
}

.fixedBinRecSampler <-
    function(buf, bin, n, tot_n)
    ## sample c('buf', 'bin') records in proportion to their
    ## reprentation, returning a list() of raw() sampled records + new
    ## buf
{
    bin <- c(buf, bin)
    env <- new.env(parent=emptyenv())
    env[["rec_n"]] <- rec_n <- .Call(.sampler_rec_counter, bin)
    if (1L >= rec_n) {
        env[["sampled"]] <- list(bin)
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
    env[["sampled"]] <-
        .Call(.sampler_rec_parser, bin, c(samp, rec_n + 1L))
    env
}

.Sampler <- setRefClass("Sampler",
    fields = list(
      con = "ANY",
      reader = "function", readerBlockSize = "integer",
      recParser = "function",
      n = "integer", saved_n = "integer", tot_n = "integer",
      records = "list", buf="raw",
      verbose = "logical"),
    methods = list(
      .sampler_msg = function(txt) {
          "display 'txt' with status information as a message()"
          s <- status()
          message(txt, " ", paste(names(s), s, sep="=", collapse=" "))
      },
      .add = function(bin, flush=FALSE) {
          ".add (incomplete) 'bin'ary stream, possibly .flush'ing buffer"
          if (verbose) .sampler_msg("sampler$.add")
          res <- recParser(buf, bin, n, tot_n)
          samp <- res[["sampled"]]
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
      .getCurrent = function() {
          "current sample _without_ .flush; incomplete and corrupt final record"
          records
      },
      .flush = function() {
          "append remaining buffer to output records"
          if (verbose) .sampler_msg("sampler$.flush")
          if (0 != length(buf)) .add(raw(), TRUE)
          .self
      },
      reset = function() {
          "reopen the connection and reset sample"
          if (verbose) .sampler_msg("sampler$reset")
          if (isOpen(con)) {
              if (verbose) .sampler_msg("sampler$reset re-open")
              s <- summary(con)
              class <- s$class
              desc <- s$description
              close(con)
              con <<- do.call(s$class, list(desc, "rb"))
          } else {
              open(con, "rb")
          }
          buf <<- raw()
          records <<- list()
          saved_n <<- tot_n <<- 0L
          .self
      },
      sample = function() {
          "read and sample all records in a connection"
          if (verbose) .sampler_msg("sampler$sample")
          reset()
          while (0 != length(bin <- reader(con, readerBlockSize)))
              .add(bin)
          .flush()
          if (verbose) .sampler_msg("sampler$sample finished")
          .self
      },
      get = function(reset=FALSE) {
          "current sample after .flush'ing; optionally invoke .reset()"
          .flush()
          result <- .getCurrent()
          if (reset) reset()
          result
      },
      status = function() {
          "report status of sampler"
          c(n=n, saved_n=saved_n, tot_n=tot_n, .buffer_len=length(buf))
      }))

##
## FastqSampler
## 

.FastqSampler <- setRefClass("FastqSampler", contains="Sampler")

FastqSampler <-
    function(con, n = 1e6, readerBlockSize=1e8, verbose=FALSE)
{
    if (length(n) != 1 || !is.finite(n) || n < 0)
        stop("'n' must be finite and >= 0")
    if (is.character(con))
        con <- file(con)
    f <- .FastqSampler$new(con=con, n=n,
                           reader=.binReader,
                           readerBlockSize=readerBlockSize,
                           recParser=.fixedBinRecSampler,
                           verbose=verbose)
}

setGeneric("yield", function(x, ...) standardGeneric("yield"))

setMethod(yield, "FastqSampler",
    function(x, ...)
{
    x$sample()
    elts <- .Call(.sampler_as_fastq, x$get())
    ShortReadQ(elts[["sread"]], elts[["quality"]], elts[["id"]], ...)
})

setMethod(show, "FastqSampler",
    function(object)
{
    cat("class:", class(object), "\n")
    cat("file:", basename(summary(object$con)$description), "\n")
    s <- object$status()
    cat("status:", paste(names(s), s, sep="=", collapse=" "), "\n")
})
