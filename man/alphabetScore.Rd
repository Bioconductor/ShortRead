\name{alphabetScore}
\alias{alphabetScore}

\title{Efficiently calculate the sum of quality scores across bases}
\description{

  This generic takes a \code{\linkS4class{QualityScore}} or
  \code{PhredQuality} object and calculates, for each read, the sum of
  the encoded nucleotide probabilities.

}
\usage{
alphabetScore(object, ...)
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{QualityScore}}.}
  \item{\dots}{Additional arguments, currently unused.}
}

\value{
  A vector of numeric values of length equal to the length of
  \code{object}.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\keyword{manip}
