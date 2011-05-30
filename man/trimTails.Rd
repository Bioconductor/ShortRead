\name{trimTails}
\alias{trimTails}
\alias{trimTails,BStringSet-method}

\title{Trim tails of reads based on nucleotide qualities}

\description{

  Short read quality typically deteriorates as cycle number
  progresses. This generic removes trailing nucleotides whose quality
  falls below a specified threshold.

}
\usage{
trimTails(object, ..., ranges=FALSE)

\S4method{trimTails}{BStringSet}(object, k, a, successive=FALSE, ..., ranges=FALSE)

}

\arguments{

  \item{object}{An object for which \code{trimTails} methods exist
    (e.g., \code{\linkS4class{ShortReadQ}} and derived classes); see
    below to discover these methods.}

  \item{k}{\code{integer(1)} describing the number of failing
    letters required to trigger trimming.}

  \item{a}{\code{character(1)} with \code{nchar(a) == 1L} giving the
    letter at or below which a nucleotide is marked as failing.}

  \item{successive}{\code{logical(1)} indicating whether failures can
    occur anywhere in the sequence, or must be successive. If
    \code{successive=FALSE}, then the k'th failed letter and subsequent
    are removed. If \code{successive=TRUE}, the first succession of k
    failed and subsequent letters are removed.}

  \item{\dots}{Additional arguments, perhaps used by methods.}

  \item{ranges}{\code{logical(1)} indicating whether the trimmed object,
    or only the ranges satisfying the trimming condition, be returned.}

}

\value{

  An instance of \code{class(object)} trimmed to contain only those
  nucleotides satisfying the trim criterion or, if \code{ranges=TRUE} an
  \code{IRanges} instance defining the ranges that would trim
  \code{object}.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
showMethods(trimTails)
}

\keyword{manip}

