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

  These generics operate on XStringSet and derived objects, and in
  particular on the short reads and their fastq-encoded quality
  scores. \code{srsort} and \code{srorder} differ from the default
  functions \code{sort} and \code{order} in that sorting is based on
  internally-defined sort order rather than, e.g., the sort order
  implied by \code{LC_COLLATE}.

}
\usage{
srorder(x, ...)
srsort(x, ...)
srduplicated(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{The object to be sorted, ordered, or duplicates identified;
    see the examples below for objects for which methods are defined.}

  \item{\dots}{Additional arguments available for use by methods;
    usually ignored.}
}
\details{

  Unlike \code{sort} and friends, the implementation does not preserve
  order of duplicated elements. Like \code{duplicated}, one element in a
  set of duplicates is marked as \code{FALSE}.

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

  \item{srorder}{An integer vector of \code{length(x)},
    identifying the indicies that wlll bring \code{x} into sorted
    order.}

  \item{srsort}{An instance of \code{x} sorted in increasing order; as
    noted.}

  \item{srduplicated}{A logical vector of \code{length(x)} indicating
    whether the indexed element is already present. Note that, like
    \code{duplicated}, subsetting \code{x} using the result returned by
    \code{!srduplicated(x)} results in one copy of all elements,
    including one copy of elements that are represented in \code{x}
    multiple times.}

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
