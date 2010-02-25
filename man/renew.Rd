\name{renewable}

\alias{renewable}
\alias{renewable,missing-method}
\alias{renewable,character-method}
\alias{renewable,.ShortReadBase-method}

\alias{renew}
\alias{renew,.ShortReadBase-method}

\title{Renew (update) a ShortRead object with new values}

\description{

  Use \code{renew} to update an object defined in \pkg{ShortRead} with
  new values. Discover update-able classes and values with
  \code{renewable}.

}
\usage{

renewable(x, \dots)
renew(x, \dots)

}

\arguments{

  \item{x}{For \code{renewable}: \code{missing}, \code{character(1)}, or
	a class defined in the \pkg{ShortRead} package. For \code{renew}:
	an instance of a class defined in the \pkg{ShortRead} package.}

  \item{\dots}{For \code{renewable}, ignored. For \code{renew}, named
	arguments identifying which parts of \code{x} are to be renewed.}

}

\details{
  
  When invoked with no arguments \code{renewable} returns a character
  vector naming classes that can be renewed.

  When invoked with a \code{character(1)} or an instance of a
  \pkg{ShortRead} class, a list of the names and values of the elements
  that can be renewed. When \code{x} is a character vector naming a
  virtual class, then each element of the returned list is a non-virtual
  descendant of that class that can be used in renewal. This is not
  fully recursive.

  \code{renew} is always invoked with the \code{x} argument being an
  instance of a class identified by \code{renewable()}. Remaining
  arguments are name-value pairs identifying the components of \code{x}
  that are to be renewed (updated). The name-value pairs must be
  consistent with \code{renewable(x)}. The resulting object is checked
  for validity. Multiple components of the object can be updated in a
  single call to \code{renew}, allowing comparatively efficient complex
  transformations.

}

\value{
  
  \code{renewable()} returns a character vector of renewable classes.

  \code{renewable(x)} returns a named list. The names correspond to
  renewable classes, and the elements of the list correspond to
  renewable components of the class.

  \code{renew(x, \dots)} returns an object of the same class as
  \code{x}, but with components of \code{x} replaced by the named values
  of \code{\dots}.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
## discovery
renewable()
renewable("AlignedRead")
renewable("QualityScore") ## instantiable classes

## example data
sp <- SolexaPath(system.file("extdata", package="ShortRead"))
ap <- analysisPath(sp)
filt <- chromosomeFilter("chr[[:digit:]+].fa")
aln <- readAligned(ap, "s_2_export.txt", "SolexaExport",
                   filter=filt)

## renew chromosomes from 'chr1.fa' to 'chr1', etc
labels <- sub("\\\\.fa", "", levels(chromosome(aln)))
renew(aln, chromosome=factor(chromosome(aln), labels=labels))

## multiple changes -- update chromosome, offset position
renew(aln, chromosome=factor(chromosome(aln), labels=labels),
      position=1L+position(aln))

## oops! invalid instances cannot be constructed
try(renew(aln, position=1:10))

}
\keyword{manip}
