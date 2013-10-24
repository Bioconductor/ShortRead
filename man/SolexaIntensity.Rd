\name{SolexaIntensity}
\alias{SolexaIntensity}
\alias{SolexaIntensityInfo}

\title{(Legacy) Construct objects of class "SolexaIntensity" and
"SolexaIntensityInfo"}

\description{

  These function constructs objects of
  \code{\linkS4class{SolexaIntensity}} and
  \code{\linkS4class{SolexaIntensityInfo}}. It will often be more
  convenient to create these objects using parsers such as
  \code{\link{readIntensities}}.

}
\usage{
SolexaIntensity(intensity=array(0, c(0, 0, 0)),
                measurementError=array(0, c(0, 0, 0)),
                readInfo=SolexaIntensityInfo(
                  lane=integer(nrow(intensity))),
                ...)
SolexaIntensityInfo(lane=integer(0),
                    tile=integer(0)[seq_along(lane)],
                    x=integer(0)[seq_along(lane)],
                    y=integer(0)[seq_along(lane)])
}

\arguments{

  \item{intensity}{A matrix of image intensity values. Successive
    columns correspond to nucleotides A, C, G, T; four successive
    columns correspond to each cycle. Typically, derived from
    \code{"_int.txt"} files.}

  \item{measurementError}{As \code{intensity}, but measuring standard
    error. Usually derived from \code{"_nse.txt"} files.}

  \item{readInfo}{An object of class \code{AnnotatedDataFrame},
    containing information described by \code{SolexaIntensityInfo}.}

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

  An object of class \code{\linkS4class{SolexaIntensity}}, or
  \code{SolexaIntensityInfo}.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{\code{\linkS4class{SolexaIntensity}}.}
\keyword{manip}
