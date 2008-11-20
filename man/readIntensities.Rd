\name{readIntensities}

\alias{readIntensities}
\alias{readIntensities,character-method}

\title{Read Solexa 'int' and 'nse' files}

\description{

  \code{readIntensities} reads \code{_int.txt} and (optionally)
  \code{_nse.txt} files in a directory into a single object. Intensities
  and their standard errors are stored as matricies; information about
  reads (lane, tile, x, and y coordinate) are recorded in an
  \code{AnnotatedDataFrame}.

}
\usage{
readIntensities(dirPath, pattern=character(0), ...,
                intExtension="_int.txt",
                nseExtension="_nse.txt",
                withNse=TRUE, verbose=FALSE)
}

\arguments{

  \item{dirPath}{Directory path or other object (e.g.,
    \code{\linkS4class{SolexaPath}} for which methods are defined.}

  \item{pattern}{A length 1 character vector representing a regular
    expression to be combined (using, e.g., \code{paste(pattern,
    intExtension, sep="")} with \code{intExtension} or
    \code{nseExtension} to match files to be summarized.}

  \item{\dots}{Additional arguments, e.g., to \code{\link{srapply}},
    used during evaluation.}

  \item{intExtension}{A length 1 character vector appended to
    \code{pattern} to select intensity files.}

  \item{nseExtension}{A length 1 character vector appended to
    \code{pattern} to select standard error (nse) files.}

  \item{withNse}{A length 1 logical vector indicating whether 'nse'
    values are to be parsed.}

  \item{verbose}{A length 1 logical vector controling information,
    useful for debuging, presented during parsing.}

}

\value{
  An object of class \code{\linkS4class{ShortReadIntensity}}.
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
fl <- system.file("extdata", package="ShortRead")
sp <- SolexaPath(fl)
int <- readIntensities(sp)
int
intensity(int)[1:5,1:5]
}

\keyword{manip}
