\name{filterFastq}

\alias{filterFastq}

\title{Filter fastq from one file to another}

\description{

  \code{filterFastq} filters reads from source to destination file(s)
  applying a filter to reads in each file. The filter can be a function
  or FilterRules instance; operations are done in a memory-efficient
  manner.

}

\usage{
filterFastq(files, destinations, ..., filter = FilterRules(),
    compress=TRUE, yieldSize = 1000000L)
}

\arguments{

  \item{files}{a character vector of valid file paths.}

  \item{destinations}{a character vector of destinations, recycled to be
    the same length as \code{files}. \code{destinations} must not
    already exist.}

  \item{...}{Additional arguments, perhaps used by a \code{filter}
    function.}

  \item{filter}{A simple function taking as it's first argument a
    \code{ShortReadQ} instance and returning a modified
    \code{ShortReadQ} instance (e.g., with records or nucleotides
    removed), or a \code{FilterRules} instance specifying which records
    are to be removed.}

  \item{compress}{A logical(1) indicating whether the file should be
    gz-compressed. The default is \code{TRUE}.}

  \item{yieldSize}{Number of fastq records processed in each call to
    \code{filter}; increase this for (marginally) more efficient I/O at
    the expense of increased memory use.}

}

\author{Martin Morgan \url{mtmorgan@fhcrc.org}}

\examples{
## path to a convenient fastq file
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")

## filter reads to keep those with GC < 0.7
fun <- function(x) {
    gc <- alphabetFrequency(sread(x), baseOnly=TRUE)[,c("G", "C")]
    x[rowSums(gc) / width(x) < .7]    
}
filterFastq(fl, tempfile(), filter=fun)

## trimEnds,character-method uses filterFastq internally
trimEnds(fl, "V", destinations=tempfile())

}
