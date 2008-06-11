\name{tables}
\alias{tables}
\alias{tables,XStringSet-method}

\title{Summarize XStringSet read frequencies}

\description{
  This generic summarizes the number of times each sequence occurs in an
  \code{\linkS4class{XStringSet}} instance.
}

\usage{
tables(x, n=50, ...)
}

\arguments{

  \item{x}{An object for which a \code{tables} method is defined.}

  \item{n}{An \code{integer(1)} value determining how many named
    sequences will be present in the \code{top} portion of the return value.}

  \item{\dots}{Additional arguments available to methods}
}

\details{

  Methods of this generic summarize the frequency with which each read
  occurs, There are two components to the summary. The reads are
  reported from most common to least common; typically a method
  parameter controls how many reads to report. Methods also return a
  pair of vectors describing how many reads were represented 1, 2,
  ... times. 

  The following methods are defined, in addition to methods described in
  class-specific documentation:
  \describe{

    \item{tables}{\code{signature(x= "XStringSet", n = 50)}:
      Apply \code{tables} to the \code{XStringSet} \code{x}.}

  }
}

\value{
  A list of length two.
  \item{top}{A named integer vector. Names correspond to sequences.
    Values are the number of times the corresponding sequence occurs in
    the \code{XStringSet}. The vector is sorted in decreasing order;
    methods typically include a parameter specifying the number of
    sequences to return.}

  \item{distribution}{a \code{data.frame} with two
    columns. \code{nOccurrences} is the number of times any particular
    sequence is represented in the set (1, 2, ...). \code{nReads} is the
    number of reads with the corresponding occurrence.}
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
showMethods("tables")
sp <- SolexaPath(system.file("extdata", package="ShortRead"))
aln <- readAligned(sp)
tables(sread(aln), n=6)
xyplot(log10(nReads)~log10(nOccurrences),
       tables(sread(aln))$distribution)
}
\keyword{manip}
