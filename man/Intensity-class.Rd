\name{Intensity-class}
\docType{class}
\alias{Intensity-class}
\alias{IntensityMeasure-class}
\alias{IntensityInfo-class}
\alias{ArrayIntensity-class}
\alias{ArrayIntensity}
% accessors
\alias{readInfo}
\alias{intensity}
\alias{measurementError}
% subset
\alias{[,IntensityMeasure,ANY,ANY,ANY-method}
\alias{[,IntensityMeasure,ANY,missing,ANY-method}
\alias{[,IntensityMeasure,missing,ANY,ANY-method}
\alias{[,IntensityMeasure,missing,missing,ANY-method}
\alias{[[,ArrayIntensity,ANY,ANY-method}
% etc.
\alias{dim,Intensity-method}
\alias{show,Intensity-method}
\alias{show,IntensityMeasure-method}

\title{"Intensity", "IntensityInfo", and "IntensityMeasure" base
  classes for short read image intensities}

\description{

  The \code{Intensity}, \code{IntensityMeasure}, and \code{IntensityInfo}
  classes represent and manipulate image intensity measures. Instances
  from the class may also contain information about measurement errors,
  and additional information about the reads from which the intensities
  are derived.

  \code{Intensity}, and \code{IntensityMeasure}, are virtual classes,
  and cannot be created directly. Classes derived from
  \code{IntensityMeasure} (e.g., \code{ArrayIntensity}) and
  \code{Intensity} (e.g., \code{\linkS4class{SolexaIntensity}}) are used
  to represent specific technologies.

}

\section{Objects from the Class}{

  \code{ArrayIntensity} objects can be created with calls of the form
  \code{ArrayIntensity(array(0, c(1,2,3)))}.

  Objects of derived classes can be created from calls such as the
  \code{\link{SolexaIntensity}} constructor, or more typically by
  parsing appropriate files (e.g., \code{\link{readIntensities}}).

}
\section{Slots}{

  Class \code{Intensity} has slots:

  \describe{

    \item{\code{readInfo}:}{Object of class \code{"IntensityInfo"}
      containing columns for the lane, tile, x, and y coordinates of the
      read.}

    \item{\code{intensity}:}{Object of class
    \code{"IntensityMeasure"} containing image intensity data for each
    read and cycle.}

    \item{\code{measurementError}:}{Object of class
      \code{"IntensityMeasure"} containing measures of image intensity
      uncertainty for each read and cycle.}

    \item{\code{.hasMeasurementError}:}{Length 1 logical variable
      indicating whether intensity standard errors are included
      (internal use only).}

  }

  Classes \code{IntensityInfo} and \code{IntensityMeasure} are virtual
  classes, and have no slots.

}
\section{Extends}{
These classes extend \code{"\linkS4class{.ShortReadBase}"}, directly.
}
\section{Methods}{

  Methods and accessor functions for \code{Intensity} include:

  \describe{

    \item{readInfo}{\code{signature(object = "Intensity")}: access the
      \code{readInfo} slot of \code{object}.}

    \item{intensity}{\code{signature(object = "Intensity")}: access the
      \code{intensity} slot of \code{object}.}

    \item{measurementError}{\code{signature(object = "Intensity")}:
      access the \code{nse} slot of \code{object}, or signal an error if
      no standard errors are available.}

    \item{dim}{\code{signature(object = "Intensity")}: return the
      dimensions (e.g., number of reads by number of cycles) represented
      by \code{object}.}

    \item{show}{\code{signature(object = "Intensity")}: provide a
      compact representation of the object.}

  }

  Subsetting \code{"["} is available for the \code{IntensityMeasure}
  class; the \code{drop} argument to \code{"["} is ignored.

  Subsetting with \code{"[["} is available for the \code{ArrayIntensity}
  class. The method accepts three arguments, corresponding to the read,
  base, and cycle(s) to be selected. The return value is the array
  (i.e., underlying data values) corresponding to the selected indices.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\link{readIntensities}}
}
\examples{
showMethods(class="Intensity", where=getNamespace("ShortRead"))
example(readIntensities)
}
\keyword{classes}
