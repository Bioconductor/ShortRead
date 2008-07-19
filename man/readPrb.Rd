\name{readPrb}

\alias{readPrb}
\alias{readPrb,character-method}

\title{Read Solexa prb files as fastq-style quality scores}

\description{

  \code{readPrb} reads all \code{_prb.txt} files in a directory into a
  single object of class \code{\linkS4class{SFastqQuality}}. Methods do
  this by identifying the maximum base call quality for each cycle and
  read, and representing this as an ASCII-encoded character string
  following Solexa conventions.

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

\value{
  An object of class \code{\linkS4class{SFastqQuality}}
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{
fl <- system.file("extdata", package="ShortRead")
sp <- SolexaPath(fl)
readPrb(sp, "s_1.*_prb.txt") # all tiles to a single file
}

\keyword{manip}
