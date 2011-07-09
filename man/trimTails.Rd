\name{trimTails}
\alias{trimTailw}
\alias{trimTailw,BStringSet-method}
\alias{trimTailw,XStringQuality-method}
\alias{trimTails}
\alias{trimTails,BStringSet-method}
\alias{trimTails,XStringQuality-method}

\title{Trim tails of reads based on nucleotide qualities}

\description{

  Short read quality typically deteriorates as cycle number
  progresses. These generics removes trailing nucleotides whose quality
  falls below a specified threshold.

}
\usage{
trimTailw(object, k, a, halfwidth, ..., ranges=FALSE)
trimTails(object, k, a, successive=FALSE, ..., ranges=FALSE)

\S4method{trimTailw}{BStringSet}(object, k, a, halfwidth, ..., alphabet, ranges=FALSE)
\S4method{trimTailw}{XStringQuality}(object, k, a, halfwidth, ..., ranges=FALSE)

\S4method{trimTails}{BStringSet}(object, k, a, successive=FALSE, ...,
    alphabet, ranges=FALSE)
\S4method{trimTails}{XStringQuality}(object, k, a, successive=FALSE, ..., ranges=FALSE)

}

\arguments{

  \item{object}{An object for which \code{trimTails} methods exist
    (e.g., \code{\linkS4class{ShortReadQ}} and derived classes); see
    below to discover these methods.}

  \item{k}{\code{integer(1)} describing the number of failing
    letters required to trigger trimming.}

  \item{a}{\code{character(1)} with \code{nchar(a) == 1L} giving the
    letter at or below which a nucleotide is marked as failing.}

  \item{halfwidth}{The half width (cycles before or after the current;
    e.g., a half-width of 5 would span 5 + 1 + 5 cycles) in which
    qualities are assessed.}

  \item{successive}{\code{logical(1)} indicating whether failures can
    occur anywhere in the sequence, or must be successive. If
    \code{successive=FALSE}, then the k'th failed letter and subsequent
    are removed. If \code{successive=TRUE}, the first succession of k
    failed and subsequent letters are removed.}

  \item{\dots}{Additional arguments, perhaps used by methods.}

  \item{alphabet}{\code{character()} (ordered low to high) letters on
    which quality scale is measured. Usually supplied internally (user
    does not need to specify). If missing, then set to ASCII characters
    0-127.}

  \item{ranges}{\code{logical(1)} indicating whether the trimmed object,
    or only the ranges satisfying the trimming condition, be returned.}

}

\details{

  \code{trimTailw} starts at the left-most nucleotide, tabulating the
  number of cycles in a window of \code{2 * halfwidth + 1} surrounding
  the current nucleotide with quality scores that fall at or below
  \code{a}. The read is trimmed at the first nucleotide for which this
  number \code{>= k}. The quality of the first or last nucleotide is
  used to represent portions of the window that extend beyond the
  sequence.

  \code{trimTails} starts at the left-most nucleotide and accumulates
  cycles for which the quality score is at or below \code{a}. The read
  is trimmed at the first location where this number \code{>= k}. With
  \code{successive=TRUE}, failing qualities must occur in strict
  succession.

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

