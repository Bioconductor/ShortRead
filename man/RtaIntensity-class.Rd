\name{RtaIntensity-class}
\docType{class}
\alias{RtaIntensity-class}

\title{Class "RtaIntensity"}
\description{

  Subclass of \code{\linkS4class{Intensity}} for representing image
  intensity data from the Illumina RTA pipeline.

}
\section{Objects from the Class}{

  Objects can be created by calls to \code{RtaIntensity} or more usually
  \code{readIntensities}.

}
\section{Slots}{

  Object of \code{RtaIntensity} have slots:

  \describe{

    \item{\code{readInfo}:}{Object of class \code{"RtaIntensityInfo"}
      representing information about each read.}

    \item{\code{intensity}:}{Object of class \code{"ArrayIntensity"}
      containing an array of intensities with dimensions read, base, and
      cycle. Nucleotide are A, C, G, T for each cycle.}

    \item{\code{measurementError}:}{Object of class
      \code{"ArrayIntensity"} containing measurement errors for each
      read, cycle, and base, with dimensions like that for
      \code{intensity}. }

    \item{\code{.hasMeasurementError}:}{Object of class
      \code{"ScalarLogical"} used internally to indicate whether
      measurement error information is included.}
  }

}
\section{Extends}{

  Class \code{"\linkS4class{SolexaIntensity}"}, directly.

  Class \code{"\linkS4class{Intensity}"}, by class "SolexaIntensity", distance 2.

  Class \code{"\linkS4class{.ShortReadBase}"}, by class "SolexaIntensity", distance 3.
}
\section{Methods}{

  Class "RtaIntensity" inherits accessor, subsetting, and display
  methods from class \code{\linkS4class{SolexaIntensity}}.

}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{\code{\linkS4class{SolexaIntensity}}, \code{\link{readIntensities}}}

\examples{
showClass("RtaIntensity")
showMethods(class="RtaIntensity", where=getNamespace("ShortRead"))
}
\keyword{classes}
