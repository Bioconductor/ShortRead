\name{ShortReadIntensity-class}
\docType{class}
\alias{ShortReadIntensity-class}
\alias{readInfo}
\alias{intensity}
\alias{nse}
\alias{show,ShortReadIntensity-method}

\title{"ShortReadIntensity" class for short read image intensities}

\description{

  This class represents and manipulates rectangular (i.e., all reads
  with the same number of nucleotides) image intensities, and
  information about their flow cell location (lane, tile, x, y
  coordinates). The class may contain measures of standard errors.

}

\section{Objects from the Class}{

  Objects of this class can be created from a call to the
  \code{\link{ShortReadIntensity}} constructor, or more typically by
  parsing appropriate files (e.g., \code{\link{readIntensities}}).

}
\section{Slots}{
  \describe{

    \item{\code{readInfo}:}{Object of class \code{"AnnotatedDataFrame"}
      containing columns for the lane, tile, x, and y coordinates of the
      read.}

    \item{\code{intensity}:}{Object of class \code{"matrix"} containing
    image intensities of each nucleotide (A, C, G, T) for each cycle
    (columns), for each read (rows).}

    \item{\code{nse}:}{Object of class \code{"matrix"} containing image
    intensity standard errors of each nucleotide (A, C, G, T) for each
    cycle (columns), for each read (rows).}

    \item{\code{.hasStandardErrors}:}{Length 1 logical variable
      indicating whether intensity standard errors are included
      (internal use only).}

  }
}
\section{Extends}{
Class \code{"\linkS4class{.ShortReadBase}"}, directly.
}
\section{Methods}{

  Methods and accessor functions for this class include:

  \describe{

    \item{readInfo}{\code{signature(object = "ShortReadIntensity")}:
      access the \code{readInfo} slot of \code{object}.}

    \item{intensity}{\code{signature(object = "ShortReadIntensity")}:
      access the \code{intensity} slot of \code{object}.}

    \item{nse}{\code{signature(object = "ShortReadIntensity")}: access
      the \code{nse} slot of \code{object}, or signal an error if no
      standard errors are available.}

    \item{show}{\code{signature(object = "ShortReadIntensity")}: provide
      a compact representation of the object.}
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\link{readIntensities}}, \code{link{ShortReadIntensity}}
}
\examples{
showMethods(class="ShortReadIntensity", where=getNamespace("ShortRead"))
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
readIntensities(sp)
}
\keyword{classes}
