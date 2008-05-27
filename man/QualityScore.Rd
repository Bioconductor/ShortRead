\name{QualityScore}
\alias{NumericQuality}
\alias{MatrixQuality}
\alias{FastqQuality}
\alias{SFastqQuality}

\title{Construct objects indicating read or alignment quality}
\description{
  Use these functions to construct quality indicators for reads or
  alignments. See \code{\linkS4class{QualityScore}} for details of
  object content and methods available for manipulating them.
}
\usage{
NumericQuality(quality = numeric(0))
MatrixQuality(quality = new("matrix"))
FastqQuality(quality=BStringSet(character(0)))
SFastqQuality(quality=BStringSet(character(0)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{quality}{An object used to initialize the data
    structure. Appropriate objects are indicated in the constructors
    above.}
}
\value{
  Constructors return objects of the corresponding class derived from
  \code{\linkS4class{QualityScore}}.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{

  \code{\linkS4class{QualityScore}}, \code{\link{readFastq}},
  \code{\link{readAligned}}

}

\examples{
nq <- NumericQuality(rnorm(20))
nq
quality(nq)
quality(nq[10:1])
}
\keyword{ manip }
