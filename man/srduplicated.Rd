\name{srduplicated}
\alias{srorder}
\alias{srrank}
\alias{srsort}
\alias{srduplicated}
% XStringSet-methods
\alias{srorder,XStringSet-method}
\alias{srrank,XStringSet-method}
\alias{srsort,XStringSet-method}
\alias{srduplicated,XStringSet-method}
\title{Order, sort, and find duplicates in XStringSet objects}
\description{

  These generics order, rank, sort, and find duplicates in short read
  objects, including fastq-encoded qualities. \code{srorder},
  \code{srrank} and \code{srsort} differ from the default functions
  \code{rank}, \code{order} and \code{sort} in that sorting is based on
  an internally-defined order rather than, e.g., the order implied by
  \code{LC_COLLATE}.

}
\usage{
srorder(x, ...)
srrank(x, ...)
srsort(x, ...)
srduplicated(x, ...)
}
\arguments{

  \item{x}{The object to be sorted, ranked, ordered, or to have
    duplicates identified; see the examples below for objects for which
    methods are defined.}

  \item{\dots}{Additional arguments available for use by methods;
    usually ignored.}
}
\details{

  Unlike \code{sort} and friends, the implementation does not preserve
  order of duplicated elements. Like \code{duplicated}, one element in
  each set of duplicates is marked as \code{FALSE}.

  \code{srrank} settles ties using the \dQuote{min} criterion described
  in \code{\link{rank}}, i.e., identical elements are ranked equal to the
  rank of the first occurence of the sorted element.

  The following methods are defined, in addition to methods described in
  class-specific documentation:
  \describe{

    \item{srsort}{\code{signature(x = "XStringSet")}:}
    \item{srorder}{\code{signature(x = "XStringSet")}:}
    \item{srduplicated}{\code{signature(x = "XStringSet")}:

      Apply \code{srorder}, \code{srrank}, \code{srsort},
      \code{srduplicated} to
      \code{\link[Biostrings:XStringSet-class]{XStringSet}} objects such
      as those returned by \code{\link{sread}}.}


    \item{srsort}{\code{signature(x = "ShortRead")}:}
    \item{srorder}{\code{signature(x = "ShortRead")}:}
    \item{srduplicated}{\code{signature(x = "ShortRead")}:

      Apply \code{srorder}, \code{srrank}, \code{srsort},
      \code{srduplicated} to
      \code{\link[Biostrings:XStringSet-class]{XStringSet}} objects to
      the \code{sread} component of \code{\linkS4class{ShortRead}} and
      derived objects.}

  }
}

\value{
  The functions return the following values:

  \item{srorder}{An integer vector the same length as \code{x},
    containing the indices that will bring \code{x} into sorted order.}

  \item{srrank}{An integer vector the same length as \code{x},
    containing the rank of each seqeunce when sorted.}

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
