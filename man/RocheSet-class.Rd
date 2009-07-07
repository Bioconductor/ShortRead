\name{RocheSet-class}
\docType{class}
\alias{RocheSet-class}
\alias{RocheSet}

\title{Roche (454) experiment-wide data container}
\description{

  This class is meant to coordinate all data in a Roche (454)
  experiment. See \code{\linkS4class{SRSet}} for additional details.

}
\section{Objects from the Class}{
  Create objects from this class using one of the \code{RocheSet}
  methods documented below
}

\section{Slots}{
  \describe{

    \item{\code{sourcePath}:}{Object of class \code{"RochePath"} The
      file system location of the data used in this experiment.}

    \item{\code{readIndex}:}{Object of class \code{"integer"} indexing
      reads included in the experiment; see \code{\linkS4class{SRSet}}
      for details on data representation in this class.}

    \item{\code{readCount}:}{Object of class \code{"integer"}
      containing the number of reads associated with each sample; see
      \code{\linkS4class{SRSet}} for details on data representation in
      this class.}

    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}
      with as many rows as there are samples, containing information on
      experimental design.}

    \item{\code{readData}:}{Object of class \code{"AnnotatedDataFrame"}
      containing as many rows as there are reads, containing information
      on each read in the experiment.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SRSet}"}, directly.
Class \code{"\linkS4class{.Roche}"}, directly.
Class \code{"\linkS4class{.ShortReadBase}"}, by class "SRSet", distance 2.
Class \code{"\linkS4class{.ShortReadBase}"}, by class ".Roche", distance 2.
}
\section{Methods}{

  No methods defined with class "RocheSet" in the signature; see
  \code{\linkS4class{SRSet}} for inherited methods.

}

\author{Michael Lawrence <mflawrence@fhcrc.org>}

\seealso{
  \code{\linkS4class{SRSet}}
}
\examples{
showClass("RocheSet")
}
\keyword{classes}
