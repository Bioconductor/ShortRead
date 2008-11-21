\name{readIntensities}

\alias{readIntensities}
\alias{readIntensities,character-method}

\title{Read Solexa 'int' and 'nse' files}

\description{

  \code{readIntensities} reads image `intensity' files (such as Solexa's
  \code{_int.txt} and (optionally) \code{_nse.txt}) in a directory into
  a single object.

}
\usage{
readIntensities(dirPath, pattern=character(0), ...)
}

\arguments{

  \item{dirPath}{Directory path or other object (e.g.,
    \code{\linkS4class{SolexaPath}} for which methods are defined.}

  \item{pattern}{A length 1 character vector representing a regular
    expression to be combined (using, e.g., \code{paste(pattern,
    intExtension, sep="")} with \code{intExtension} or
    \code{nseExtension} to match files to be summarized.}

  \item{\dots}{Additional arguments used by methods.}

}

\details{

  Additional methods are defined on specific classes, see, e.g.,
  \code{\linkS4class{SolexaPath}}.

  The \code{readIntensities,character-method} contains an arugment
  \code{type} that determines how intensities are parsed. Use the
  \code{type} argument to \code{readIntensities,character-method}, as
  follows:

  \itemize{

    \item{\code{type="SolexaIntensity"}}{

      Intensities are read from Solexa \code{_int.txt} and
      \code{_nse.txt}-style files; see \code{\linkS4class{SolexaPath}}
      for details.

    }
  }
}

\value{
  An object derived from class \code{\linkS4class{Intensity}}.
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
fl <- system.file("extdata", package="ShortRead")
sp <- SolexaPath(fl)
int <- readIntensities(sp)
int
intensity(int)[1,,] 		# one read
head(rowMeans(intensity(int)))	# treat as 'matrix'
head(pData(readInfo(int)))
}

\keyword{manip}
