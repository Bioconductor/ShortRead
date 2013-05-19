\name{trimTails}
\alias{trimTailw}
\alias{trimTailw,character-method}
\alias{trimTailw,BStringSet-method}
\alias{trimTailw,XStringQuality-method}

\alias{trimTails}
\alias{trimTails,character-method}
\alias{trimTails,BStringSet-method}
\alias{trimTails,XStringQuality-method}

\alias{trimEnds}
\alias{trimEnds,character-method}
\alias{trimEnds,XStringSet-method}
\alias{trimEnds,XStringQuality-method}
\alias{trimEnds,FastqQuality-method}
\alias{trimEnds,ShortRead-method}
\alias{trimEnds,ShortReadQ-method}

\title{Trim ends of reads based on nucleotides or qualities}

\description{

  These generic functions remove leading or trailing nucleotides or
  qualities. \code{trimTails} and \code{trimTailw} remove low-quality
  reads from the right end using a sliding window (\code{trimTailw}) or
  a tally of (successive) nucleotides falling at or below a quality
  threshold (\code{trimTails}). \code{trimEnds} takes an alphabet of
  characters to remove from either left or right end.

}
\usage{
## S4 methods for 'ShortReadQ', 'FastqQuality', or 'SFastqQuality'
trimTailw(object, k, a, halfwidth, ..., ranges=FALSE)
trimTails(object, k, a, successive=FALSE, ..., ranges=FALSE)
trimEnds(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
    ..., ranges=FALSE)

\S4method{trimTailw}{BStringSet}(object, k, a, halfwidth, ..., alphabet, ranges=FALSE)
\S4method{trimTails}{BStringSet}(object, k, a, successive=FALSE, ...,
    alphabet, ranges=FALSE)

\S4method{trimTailw}{character}(object, k, a, halfwidth, ..., destinations, ranges=FALSE)
\S4method{trimTails}{character}(object, k, a, successive=FALSE, ..., destinations, ranges=FALSE)
\S4method{trimEnds}{character}(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
    ..., destinations, ranges=FALSE)
}

\arguments{

  \item{object}{An object (e.g., \code{\linkS4class{ShortReadQ}} and
    derived classes; see below to discover these methods) or character
    vector of fastq file(s) to be trimmed. }

  \item{k}{\code{integer(1)} describing the number of failing
    letters required to trigger trimming.}

  \item{a}{For \code{trimTails} and \code{trimTailw}, a
    \code{character(1)} with \code{nchar(a) == 1L} giving the letter at
    or below which a nucleotide is marked as failing.

    For \code{trimEnds} a \code{character()} with all \code{nchar() ==
    1L} giving the letter at or below which a nucleotide or quality
    scores marked for removal.}

  \item{halfwidth}{The half width (cycles before or after the current;
    e.g., a half-width of 5 would span 5 + 1 + 5 cycles) in which
    qualities are assessed.}

  \item{successive}{\code{logical(1)} indicating whether failures can
    occur anywhere in the sequence, or must be successive. If
    \code{successive=FALSE}, then the k'th failed letter and subsequent
    are removed. If \code{successive=TRUE}, the first succession of k
    failed and subsequent letters are removed.}

  \item{left, right}{\code{logical(1)} indicating whether trimming is
    from the left or right ends.}

  \item{relation}{\code{character(1)} selected from the argument values,
    i.e., \dQuote{<=} or \dQuote{==} indicating whether all letters at
    or below the \code{alphabet(object)} are to be removed, or only
    exact matches.}

  \item{\dots}{Additional arguments, perhaps used by methods.}

  \item{destinations}{For \code{object} of type \code{character()}, an
    equal-length vector of destination files. Files must not already
    exist.}

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

  \code{trimEnds} examines the \code{left}, \code{right}, or both ends
  of \code{object}, marking for removal letters that correspond to
  \code{a} and \code{relation}. The \code{trimEnds,ShortReadQ-method}
  trims based on quality.

  \code{ShortReadQ} methods operate on quality scores; use
  \code{sread()} and the \code{ranges} argument to trim based on
  nucleotide (see examples).

  \code{character} methods transform one or several fastq files to new
  fastq files, applying trim operations based on quality scores; use
  \code{filterFastq} with your own \code{filter} argument to filter on
  nucleotides.

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

sp <- SolexaPath(system.file('extdata', package='ShortRead'))
rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")

## remove leading / trailing quality scores <= 'I'
trimEnds(rfq, "I")
## remove leading / trailing 'N's
rng <- trimEnds(sread(rfq), "N", relation="==", ranges=TRUE)
narrow(rfq, start(rng), end(rng))
## remove leading / trailing 'G's or 'C's
trimEnds(rfq, c("G", "C"), relation="==")

}

\keyword{manip}

