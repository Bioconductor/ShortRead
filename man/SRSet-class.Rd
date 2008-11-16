\name{SRSet-class}
\docType{class}
\alias{SRSet-class}
\alias{readCount}
\alias{readData}
\alias{readIndex}
\alias{sourcePath}
\alias{experimentPath,SRSet-method}
\alias{show,SRSet-method}
\alias{detail,SRSet-method}

\title{A base class for Roche experiment-wide data}
\description{
  This class coordinates phenotype (sample) and sequence data, primarily
  as used on the Roche platform.

  Conceptually, this class has reads from a single experiment
  represented as a long vector, ordered by sample. The \code{readCount}
  slot indicates the number of reads in each sample, so that the sum of
  \code{readCount} is the total number of reads in the experiment. The
  \code{readIndex} field is a light-weight indicator of which reads from
  all those available that are currently referenced by the \code{SRSet}.

}

\section{Objects from the Class}{

  Objects of this class are not usually created directly, but instead
  are created by a derived class, e.g., \code{\linkS4class{RocheSet}}.

}
\section{Slots}{
  \describe{

    \item{\code{sourcePath}:}{Object of class \code{"ExperimentPath"},
      containing the directory path where sequence files can be found.}

    \item{\code{readIndex}:}{Object of class \code{"integer"} indicating
    specific sequences included in the experiment.}

    \item{\code{readCount}:}{Object of class \code{"integer"} containing
      the number of reads in each sample included in the experiment. The
      sum of this vector is the total number of reads.}

    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}
      describing each sample in the experiment. The number of rows of
      \code{phenoData} equals the number of elements in
      \code{readCount}.}

    \item{\code{readData}:}{Object of class \code{"AnnotatedDataFrame"}
      containing annotations on all reads.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{.ShortReadBase}"}, directly.
}
\section{Methods}{
  \describe{
    \item{experimentPath}{\code{signature(object = "SRSet")}: return the
      \code{\linkS4class{ExperimentPath}} associated with this object.}

    \item{readCount}{\code{signature(object="SRSet")}:}
    \item{readIndex}{\code{signature(object="SRSet")}:}
    \item{readData}{\code{signature(object="SRSet")}:}
    \item{sourcePath}{\code{signature(object="SRSet")}: Retrieve the
      corresponding slot from \code{object}.}
    

    \item{show}{\code{signature(object = "SRSet")}: display the contents
      of this object.}

    \item{detail}{\code{signature(object = "SRSet")}: provide more
      extensive information on the object.}
  }
}

\author{Michael Lawrence <mflawrence@fhcrc.org>}

\examples{
showClass("SRSet")
}
\keyword{classes}
