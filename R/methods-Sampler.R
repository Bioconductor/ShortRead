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

.FastqFile <- setRefClass("FastqFile",
    fields=list(
      con = "ANY",
      reader = "function", readerBlockSize = "integer",
      recParser = "function",
      n = "integer", saved_n = "integer", tot_n = "integer",
      records = "list", buf="raw",
      verbose = "logical"),
    methods=list(
      finalize = function() {
          ## isOpen fails after(close(con))...
          tryCatch(if (isOpen(con)) close(con),
                   error=function(...) {})
          .self
      },
      .fastqfile_msg = function(txt) {
          "display 'txt' with status information as a message()"
          s <- status()
          message(txt, " ", paste(names(s), s, sep="=", collapse=" "))
      },
      .flush = function() {
          "append remaining buffer to output records"
          if (verbose) .fastqfile_msg("FastqFile$.flush")
          if (0 != length(buf)) .add(raw(), TRUE)
          .self
      },
      reset = function() {
          "reopen the connection and reset sample"
          if (verbose) .fastqfile_msg("FastqFile$reset")
          if (isOpen(con)) {
              if (verbose) .fastqfile_msg("FastqFile$reset re-open")
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
      get = function(reset=FALSE) {
          "current records after .flush'ing; optionally invoke .reset()"
          .flush()
          result <- records
          if (reset) reset()
          result
      },
      status = function() {
          "report status of FastqFile"
          c(n=n, saved_n=saved_n, tot_n=tot_n, .buffer_len=length(buf))
      },
      show = function() {
          cat("class:", class(.self), "\n")
          cat("file:", basename(summary(.self$con)$description), "\n")
          s <- .self$status()
          cat("status:", paste(names(s), s, sep="=", collapse=" "), "\n")
      }))

setGeneric("yield", function(x, ...) standardGeneric("yield"))

## 
## FastqSampler
## 

.FastqSampler <- setRefClass("FastqSampler", contains="FastqFile",
    methods = list(
      .add = function(bin, flush=FALSE) {
          ".add (incomplete) 'bin'ary stream, possibly .flush'ing buffer"
          if (verbose) .fastqfile_msg("FastqSampler$.add")
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
          if (verbose) .fastqfile_msg("FastqSampler$.yield")
          reset()
          while (0 != length(bin <- reader(con, readerBlockSize)))
              .add(bin)
          .flush()
          if (verbose) .fastqfile_msg("FastqSampler$sample finished")
          .self
      }))

FastqSampler <-
    function(con, n = 1e6, readerBlockSize=1e8, verbose=FALSE)
{
    if (length(n) != 1 || !is.finite(n) || n < 0)
        stop("'n' must be finite and >= 0")
    if (is.character(con))
        con <- file(con)
    .FastqSampler$new(con=con, n=as.integer(n),
                      reader=.binReader,
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

##
## FastqStreamer
##

.FastqStreamer <- setRefClass("FastqStreamer", contains="FastqFile",
    methods = list(
      .add = function(bin, flush=FALSE) {
          ".add (incomplete) 'bin'ary stream, possibly .flush'ing buffer"
          if (verbose) .fastqfile_msg("FastqSampler$.add")
          res <- recParser(buf, bin, Inf, tot_n)
          samp <- res[["parsed_bin"]]
          if (flush)
              buf <<- raw()
          if (length(samp))
              records <<- c(records, samp)
          saved_n <<- length(records)
          tot_n <<- tot_n + res[["rec_n"]]
          .self
      },
      get = function() {
          "at most 'n' complete records"
          if (verbose) .fastqfile_msg("FastqStreamer$get")

          while (n > length(records) &&
                 0 != length(bin <- reader(con, readerBlockSize)))
              .add(bin)
          if (n > length(records))
              .flush()                  # last record
          if (n < length(records)) {
              idx <- 1:n
              res <- records[idx]
              records <<- records[-idx]
          } else {
              res <- records
              records <<- list()
          }
          res
      }))

FastqStreamer <-
    function(con, n=1e6, readerBlockSize=1e8, verbose=FALSE)
{
    if (length(n) != 1 || !is.finite(n) || n < 0)
        stop("'n' must be finite and >= 0")
    if (is.character(con))
        con <- file(con, "rb")
    .FastqStreamer$new(con=con, n=as.integer(n), tot_n=0L, saved_n=0L,
                       reader=.binReader,
                       readerBlockSize=as.integer(readerBlockSize),
                       recParser=.fixedBinRecSampler,
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
