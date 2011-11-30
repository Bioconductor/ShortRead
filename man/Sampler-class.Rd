\name{FastqFile-class}
\docType{class}
% Classes
\alias{ShortReadFile-class}
\alias{FastqFile-class}
\alias{FastqFileReader-class}
\alias{FastqSampler-class}
\alias{FastqStreamer-class}
\alias{FastqFileList-class}

% ShortReadFile / FastqFile / FastqFileList
\alias{FastqFile}
\alias{FastqFileList}
\alias{open.ShortReadFile}
\alias{close.ShortReadFile}
\alias{readFastq,FastqFile-method}

% FastqFileReader
\alias{yield,FastqFileReader-method}

% FastqSampler
\alias{FastqSampler}
\alias{yield}
\alias{yield,FastqSampler-method}

% FastqStreamer
\alias{FastqStreamer}
\alias{yield,FastqStreamer-method}

\title{Sampling and streaming records from fastq files}

\description{

  \code{FastqFile} represents a path and connection to a fastq
  file. \code{FastqFileList} is a list of such connections.

  \code{FastqSampler} draws a subsample from a fastq file. \code{yield}
  is the method used to extract the sample from the \code{FastqSampler}
  instance; a short illustration is in the example below.

  \code{FastqStreamer} draws successive subsets from a fastq
  file. Iterating over the stream will eventually yield the entire
  file; a short illustration is in the example below.

}

\usage{
## FastqFile and FastqFileList
FastqFile(con, ...)
FastqFileList(...)
\S3method{open}{ShortReadFile}(con, ...)
\S3method{close}{ShortReadFile}(con, ...)
\S4method{readFastq}{FastqFile}(dirPath, pattern=character(), ...)

## FastqSampler and FastqStreamer
FastqSampler(con, n=1e6, readerBlockSize=1e8, verbose=FALSE)
FastqStreamer(con, n=1e6, readerBlockSize=1e8, verbose=FALSE)
yield(x, ...)
}

\arguments{

  \item{con, dirPath}{A character string naming a connection, or (for
    \code{con}) an R connection (e.g., \code{file}, \code{gzfile}).}

  \item{n}{The size of the sample (number of records) to be drawn.}

  \item{readerBlockSize}{The number of bytes or characters to be read at
    one time; smaller \code{readerBlockSize} reduces memory requirements
    but is less efficient.}

  \item{verbose}{Display progress.}

  \item{x}{An instance from the \code{FastqSampler} or
  \code{FastqStreamer} class.}

  \item{...}{Additional arguments. For \code{FastqFileList}, this can
    either be a single character vector of paths to fastq files, or
    several instances of \code{FastqFile} objects.}

  \item{pattern}{Ignored.}

}

\section{Objects from the class}{

  Available classes include:

  \describe{

    \item{\code{FastqFile}}{A file path and connection to a fastq file.}

    \item{\code{FastqFileList}}{A list of \code{FastqFile} instances.}

    \item{\code{FastqSampler}}{Uniformly sample records from a fastq
      file.}

    \item{\code{FastqStreamer}}{Iterate over a fastq file, returning
      successive parts of the file.}

  }
}

\section{Methods}{

  The following methods are available to users:

  \describe{

    \item{\code{readFastq,FastqFile-method}:}{see also
      \code{?\link{readFastq}}.}

    \item{\code{writeFastq,ShortReadQ,FastqFile-method}:}{see also
      \code{?\link{writeFastq}}, 
      \code{?"writeFastq,ShortReadQ,FastqFile-method"}.}

    \item{\code{yield}:}{Draw a single sample from the
      instance. Operationally this requires that the underlying data
      (e.g., file) represented by the \code{Sampler} instance be
      visited; this may be time consuming.}

  }

}

\seealso{

  \code{\link{readFastq}}, \code{\link{writeFastq}}, \code{\link{yield}}.

}

\examples{
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")

f <- FastqFile(fl)
rfq <- readFastq(f)

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
    print(length(fq))
}

## Internal fields, methods, and help; for developers
ShortRead:::.FastqSampler_g$methods()
ShortRead:::.FastqSampler_g$fields()
ShortRead:::.FastqSampler_g$help("yield")

}
