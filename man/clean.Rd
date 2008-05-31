\name{clean}
\alias{clean}
\alias{clean,DNAStringSet-method}

\title{Remove sequences with ambiguous nucleotides from short read classes}

\description{

  Short reads may contain ambiguous base calls (i.e., IUPAC symbols
  different from A, T, G, C). This generic removes all sequences
  containing 1 or more ambiguous bases.

}
\usage{
clean(object, ...)
}
\arguments{

  \item{object}{An object for which \code{clean} methods exist; see
    below to discover these methods.} 

  \item{\dots}{Additional arguments, perhaps used by methods.}
}
\details{
  The following method is defined, in addition to methods described in
  class-specific documentation:
  \describe{

    \item{clean}{\code{signature(x = "DNAStringSet")}:
      Remove all sequences containin non-base (A, C, G, T) IUPAC
      symbols.}

  }
}
\value{

  An instance of \code{class(object)}, containing only sequences with
  non-redundant nucleotides.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
showMethods('clean')
}
\keyword{manip}

