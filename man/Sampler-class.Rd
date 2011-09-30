\name{FastqFile-class}
\docType{class}
\alias{FastqFile-class}

% FastqSampler
\alias{FastqSampler-class}
\alias{FastqSampler}
\alias{yield}
\alias{yield,FastqSampler-method}

% FastqStreamer
\alias{FastqStreamer-class}
\alias{FastqStreamer}
\alias{yield,FastqStreamer-method}

\title{Sampling and streaming records from fastq files}

\description{

  \code{FastqSampler} draws a subsample from a fastq file. \code{yield}
  is the method used to extract the sample from the \code{FastqSampler}
  instance; a short illustration is in the example below.

  \code{FastqStreamer} draws successive subsets from a fastq
  file. Iterating over the stream will eventually yield the entire
  file; a short illustration is in the example below.

}

\usage{
FastqSampler(con, n=1e6, readerBlockSize=1e8, verbose=FALSE)
FastqStreamer(con, n=1e6, readerBlockSize=1e8, verbose=FALSE)
yield(x, ...)
}

\arguments{

  \item{con}{A character string naming a connection, or a connection.}

  \item{n}{The size of the sample (number of records) to be drawn.}

  \item{readerBlockSize}{The number of bytes or characters to be read at
    one time; smaller \code{readerBlockSize} reduces memory requirements
    but is less efficient.}

  \item{verbose}{Display progress.}

  \item{x}{An instance from the \code{FastqSampler} or
  \code{FastqStreamer} class.}

  \item{...}{Additional arguments; currently none.}

}

\section{Objects from the class}{

  Available classes include:

  \describe{

    \item{\code{FastqSampler}}{Uniformly sample records from a fastq
      file.}

    \item{\code{FastqStreamer}}{Iterate over a fastq file, returning
      successive parts of the file.}

  }
}

\section{Methods}{

  The following methods are available to users:

  \describe{

    \item{\code{yield}:}{Draw a single sample from the
      instance. Operationally this requires that the underlying data
      (e.g., file) represented by the \code{Sampler} instance be
      visited; this may be time consuming.}

  }

}

\seealso{
  \code{\link{FastqQuality}} for instance construction.
  \code{\link{yield}} for generic and method description.
}

\examples{
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")

f <- FastqSampler(fl, 50)
yield(f)    # sample of size n=50
yield(f)    # independent sample of size 50

f <- FastqStreamer(fl, 50)
yield(f)    # records 1 to 50
yield(f)    # records 51 to 100

## iterating over an entire file
f <- FastqStreamer(fl, 50)
while (length(fq <- yield(f))) {
    ## do work here
}

## Internal fields, methods, and help; for developers
ShortRead:::.FastqSampler$methods()
ShortRead:::.FastqSampler$fields()
ShortRead:::.FastqSampler$help("get")
## Internal -- sampled records as list of raw();
##   yield() knows how to translate these to ShortReadQ
f$status()
recs <- f$get()
str(recs[1:5])
}
