\name{RtaIntensity}
\alias{RtaIntensity}

\title{Construct objects of class "RtaIntensity"}

\description{

  \code{\linkS4class{RtaIntensity}} objects contain Illumina image
  intensity measures created by the RTA pipeline. It will often be more
  convenient to create this object using \code{\link{readIntensities}}.

}
\usage{
RtaIntensity(intensity=array(0, c(0, 0, 0)),
             measurementError=array(0, c(0, 0, 0)),
             readInfo=SolexaIntensityInfo(
               lane=integer()[seq_len(nrow(intensity))]),
             ...)
}

\arguments{

  \item{intensity}{A matrix of image intensity values. Successive
    columns correspond to nucleotides A, C, G, T; four successive
    columns correspond to each cycle. Typically, derived from
    \code{"_int.txt"} files.}

  \item{measurementError}{As \code{intensity}, but measuring standard
    error. Usually derived from \code{"_nse.txt"} files.}

  \item{readInfo}{An object of class \code{AnnotatedDataFrame},
    containing information described by \code{RtaIntensityInfo}.}

  \item{...}{Additional arguments, not currently used.}
}
\value{

  An object of class \code{\linkS4class{RtaIntensity}}.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{\code{\linkS4class{RtaIntensity}},
  \code{\link{readIntensities}}.}

\examples{
rta <- RtaIntensity(array(runif(60), c(5,4,3)))
intensity(rta)
## subsetting, access, and coercion
as(intensity(rta)[1:2,,], "array")
}
\keyword{manip}
