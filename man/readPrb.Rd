\name{readPrb}

\alias{readPrb}
\alias{readPrb,character-method}

\title{Read Solexa prb files as fastq-style quality scores}

\description{

  \code{readPrb} reads all \code{_prb.txt} files in a directory into a
  single object. Most methods (see details) do this by identifying the
  maximum base call quality for each cycle and read, and representing
  this as an ASCII-encoded character string.

}
\usage{
readPrb(dirPath, pattern = character(0), ...)
}

\arguments{

  \item{dirPath}{Directory path or other object (e.g.,
    \code{\linkS4class{SolexaPath}} for which methods are defined.}

  \item{pattern}{Regular expression matching names of \code{_prb} files
    to be summarized.}

  \item{\dots}{Additional arguments, e.g., to \code{\link{srapply}},
    used during evaluation.}

}

\details{

  The \code{readPrb,character-method} contains an arugment \code{as}
  that determines the value of the returned object, as follows.

  \itemize{

    \item{\code{as="SolexaEncoding"}}{

      The ASCII encoding of the maximum per cycle and read quality score
      is encoded using Solexa conventions.

    }

    \item{\code{as="FastqEncoding"}}{

      The ASCII encoding of the maximum per cycle and read quality score
      is encoded using Fastq conventions, i.e., \code{!} has value 0.

    }

    \item{\code{as="IntgerEncoding"}}{

      The maximum per cycle and read quality score is returned as a in
      integer value. Values are collated into a matrix with number of
      rows equal to number of reads, and number of columns equal to
      number of cycles.

    }

    \item{\code{as="array"}}{

      The quality scores are \emph{not} summarized; the return value is
      an integer array with dimensions corresponding to reads,
      nucleotids, and cycles.

    }

  }
}

\value{
  An object of class \code{\linkS4class{QualityScore}}, or an integer matrix.
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
fl <- system.file("extdata", package="ShortRead")
sp <- SolexaPath(fl)
readPrb(sp, "s_1.*_prb.txt") # all tiles to a single file
}

\keyword{manip}
