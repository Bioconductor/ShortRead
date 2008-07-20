\name{Utilites}
\alias{polyn}

\title{Utilities for common, simple operations}

\description{

  These functions perform a variety of simple operations.

}
\usage{
polyn(nucleotides, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nucleotides}{A character vector with all elements having exactly
    1 character, typically from the IUPAC alphabet.}
  \item{n}{An \code{integer(1)} vector.}
}
\details{

  \code{polyn} returns a character vector with each element having
  \code{n} characters. Each element contains a single nucleotide. Thus
  \code{polyn("A", 5)} returns \code{AAAAA}.

}
\value{

  \code{polyn} returns a character vector of length \code{length(nucleotide)}

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
polyn(c("A", "N"), 35)
}

\keyword{manip}
