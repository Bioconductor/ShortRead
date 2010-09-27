\name{Sampler-class}
\docType{class}
\alias{Sampler-class}
\alias{FastqSampler-class}
\alias{FastqSampler}
\alias{yield}
\alias{yield,FastqSampler-method}
\alias{show,FastqSampler-method}

\title{Sampling records from fastq files}

\description{

  The \code{Sampler} class represents a subsample of records from a
  file. \code{FastqSampler} is an implementation to sample from a fastq
  file. \code{yield} is the method used to extract the sample from the
  \code{Sampler} or \code{FastqSampler} class; a short illustration is
  in the example below.

}

\usage{
FastqSampler(con, n=1e6, readerBlockSize=1e8, verbose=FALSE)
yield(x, ...)
\S4method{yield}{FastqSampler}(x, ...)
}

\arguments{

  \item{con}{A character string naming a connection, or a connection.}

  \item{n}{The size of the sample (number of records) to be drawn.}

  \item{readerBlockSize}{The number of bytes or characters to be read at
    one time; smaller \code{readerBlockSize} reduces memory requirements
    but is less efficient.}

  \item{verbose}{Display progress.}

  \item{x}{An instance from a class extending the \code{Sampler}-class.}

  \item{...}{Additional arguments; currently none.}

}

\section{Objects from the class}{

  Available \code{Sampler} classes include:

  \describe{

    \item{\code{Sampler}}{Base class; requires implementation.}

    \item{\code{FastqSampler}}{Uniformly sample records from a fastq
      file. See the \code{\link{FastqSampler}} constructor and help
      pages for methods mentioned below. \code{FastqSampler} extends
      \code{Sampler}. }

  }
}

\section{Methods}{

  The following methods are available to users:

  \describe{

    \item{\code{yield}:}{Draw a single sample from the
      instance. Operationally this requires that the underlying data
      (e.g., file) represented by the \code{Sampler} instance be
      visited; this may be time consuming.}

    \item{\code{show}:}{Display summary information about the instance.}

  }

  \code{Sampler} and derived classes are `reference' classes. The
  intended implementation is that users do not access the methods and
  fields of the class directly, but instead use S4 methods defined on
  the class. Nonetheless, methods and fields are available by accessing
  the class definitions. Field and method documentation are as described
  in \code{?ReferenceClasses}.

}

\seealso{
  \code{\link{FastqQuality}} for instance construction.
  \code{\link{yield}} for generic and method description.
}

\examples{
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")
f <- FastqSampler(fl, 50, verbose=TRUE)
yield(f)    # sample of size n=50
yield(f)    # independent sample of size 50


## Internal fields, methods, and help; for developers
ShortRead:::.Sampler$methods()
ShortRead:::.Sampler$fields()
ShortRead:::.Sampler$help("get")
## Internal -- sampled records as list of raw();
##   yield() knows how to translate these to ShortReadQ
f$status()
recs <- f$get()
str(recs[1:5])
}
