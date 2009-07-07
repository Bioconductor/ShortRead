\name{SolexaIntensity-class}
\docType{class}
\alias{SolexaIntensity-class}
\alias{SolexaIntensityInfo-class}
% subset
\alias{[,SolexaIntensity,ANY,ANY-method}
\alias{[,SolexaIntensity,ANY,ANY,ANY-method}
\alias{[,SolexaIntensity,ANY,missing,ANY-method}
\alias{[,SolexaIntensity,missing,ANY,ANY-method}
\alias{[,SolexaIntensity,missing,missing,ANY-method}

\title{Classes "SolexaIntensity" and "SolexaIntensityInfo"}
\description{

  Instances of \code{\linkS4class{Intensity}} and
  \code{\linkS4class{IntensityInfo}} for representing image intensity
  data from Solexa experiments.

}
\section{Objects from the Class}{

  Objects can be created by calls to \code{SolexaIntensityInfo} or
  \code{SolexaIntensity}, or more usually \code{readIntensities}.

}
\section{Slots}{

  Object of \code{SolexaIntensity} have slots:

  \describe{

    \item{\code{readInfo}:}{Object of class \code{"SolexaIntensityInfo"}
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

  Object of \code{SolexaIntensityInfo} 

  \describe{

    \item{data}{Object of class \code{"data.frame"}, inherited from
      \code{AnnotatedDataFrame}.}
    
    \item{varMetadata}{Object of class \code{"data.frame"}, inherited
      from \code{AnnotatedDataFrame}.}

    \item{dimLabels}{Object of class \code{"character"}, inherited
      from \code{AnnotatedDataFrame}.}

    \item{\code{.__classVersion__}}{Object of class \code{"Versions"},
      inherited from \code{AnnotatedDataFrame}.}

    \item{.init}{Object of class \code{"ScalarLogical"}, used internally
      to indicate whether the user initialized this object.}

  }

}
\section{Extends}{

  Class \code{SolexaIntensity}:

  Class \code{"\linkS4class{Intensity}"}, directly.
  Class \code{"\linkS4class{.ShortReadBase}"}, by class "Intensity", distance 2.

  Class \code{SolexaIntensityInfo}:

  Class \code{"\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}"}, directly
  Class \code{"\linkS4class{IntensityInfo}"}, directly
  Class \code{"\link[Biobase:class.Versioned]{Versioned}"}, by class "AnnotatedDataFrame", distance 2
  Class \code{"\linkS4class{.ShortReadBase}"}, by class "IntensityInfo", distance 2  Class \code{"\linkS4class{IntensityInfo}"}, directly.
}
\section{Methods}{

  Class "SolexaIntensity" inherits accessor and display
  methods from class \code{\linkS4class{Intensity}}. Additional methods include:

  \describe{

    \item{\code{[}}{\code{signature(x = "SolexaIntensity", i="ANY", j="ANY", k="ANY")}: 

      Selects the ith read, jth nucleotide, and kth cycle. Selection is
      coordinated across intensity, measurement error, and read
      information.}

  }

  Class "SolexaIntensityInfo" inherits accessor, subsetting, and display
  methods from class \code{\linkS4class{IntensityInfo}} and
  \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}.

}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{\code{\link{readIntensities}}}

\examples{
showClass("SolexaIntensity")
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
int <- readIntensities(sp)
int 				# SolexaIntensity
readInfo(int)			# SolexaIntensityInfo
int[1:5,,]			# read 1:5
}
\keyword{classes}
