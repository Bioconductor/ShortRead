\name{srduplicated}
\alias{srorder}
\alias{srsort}
\alias{srduplicated}
% XStringSet-methods
\alias{srorder,XStringSet-method}
\alias{srsort,XStringSet-method}
\alias{srduplicated,XStringSet-method}
\title{Order, sort, and find duplicates in XStringSet objects}
\description{

  These generics sort, order, and find duplicates in short read objects,
  including fastq-encoded qualities. \code{srsort} and \code{srorder}
  differ from the default functions \code{sort} and \code{order} in that
  sorting is based on an internally-defined order rather than, e.g., the
  order implied by \code{LC_COLLATE}.

}
\usage{
srorder(x, ...)
srsort(x, ...)
srduplicated(x, ...)
}
\arguments{

  \item{x}{The object to be sorted, ordered, or to have duplicates
    identified; see the examples below for objects for which methods are
    defined.}

  \item{\dots}{Additional arguments available for use by methods;
    usually ignored.}
}
\details{

  Unlike \code{sort} and friends, the implementation does not preserve
  order of duplicated elements. Like \code{duplicated}, one element in
  each set of duplicates is marked as \code{FALSE}.

  The following methods are defined, in addition to methods described in
  class-specific documentation:
  \describe{

    \item{srsort}{\code{signature(x = "XStringSet")}:}
    \item{srorder}{\code{signature(x = "XStringSet")}:}
    \item{srduplicated}{\code{signature(x = "XStringSet")}:

      Apply \code{srsort}, \code{srorder}, \code{srduplicated} to
      \code{\linkS4class{XStringSet}} objects such as those returned by
      \code{\link{sread}}.}

  }
}

\value{
  The functions return the following values:

  \item{srorder}{An integer vector the same length as \code{x},
    containing the indicies that wlll bring \code{x} into sorted order.}

  \item{srsort}{An instance of \code{x} in sorted order.}

  \item{srduplicated}{A logical vector the same length as \code{x}
    indicating whether the indexed element is already present. Note
    that, like \code{duplicated}, subsetting \code{x} using the result
    returned by \code{!srduplicated(x)} includes one representative from
    each set of duplicates.}

}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
showMethods("srsort")
showMethods("srorder")
showMethods("srduplicated")

sp <- SolexaPath(system.file('extdata', package='ShortRead'))
rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")

sum(srduplicated(sread(rfq)))
srsort(sread(rfq))
srsort(quality(rfq))
}
\keyword{manip}
