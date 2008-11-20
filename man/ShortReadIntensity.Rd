\name{ShortReadIntensity}
\alias{ShortReadIntensity}
\alias{ShortReadIntensityInfo}

\title{Construct objects of class "ShortReadIntensity"}
\description{

  This function constructs objects of
  \code{\linkS4class{ShortReadIntensity}}. It will often be more
  convenient to create \code{ShortReadIntensity} objects using parsers
  such as \code{\link{readIntensities}}.

}
\usage{
ShortReadIntensity(intensity=matrix(0, 0, 0),
                   nse=matrix(0, 0, 0),
                   readInfo=ShortReadIntensityInfo(
                     lane=integer(nrow(intensity))),
                   ...)
ShortReadIntensityInfo(lane=integer(0),
                       tile=integer(length(lane)),
                       x=integer(length(lane)),
                       y=integer(length(lane)))
}

\arguments{

  \item{intensity}{A matrix of image intensity values. Successive columns
    correspond to nucleotides A, C, G, T; four successive columns
    correspond to each cycle.}

  \item{nse}{As \code{intensity}, but measuring standard error.}

  \item{readInfo}{An object of class \code{AnnotatedDataFrame},
    containing information described by \code{ShortReadIntensityInfo}.}

  \item{lane}{An integer vector giving the lane from which each read is
    derived.}

  \item{tile}{An integer vector giving the tile from which each read is
    derived.}

  \item{x}{An integer vector giving the tile-local x coordinate of the
    read from which each read is derived.}

  \item{y}{An integer vector giving the tile-local y coordinate of the
    read from which each read is derived.}

  \item{...}{Additional arguments, not currently used.}
}
\value{

  An object of class \code{\linkS4class{ShortReadIntensity}}, or (for
  \code{ShortReadIntensityInfo}) an \code{AnnotatedDataFrame}.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{\code{\linkS4class{ShortReadIntensity}}.}
\keyword{manip}
