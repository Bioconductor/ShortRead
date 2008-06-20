\name{AlignedDataFrame-class}
\docType{class}
\alias{AlignedDataFrame-class}

\title{
  "AlignedDataFrame" representing alignment annotations as a data frame
}
\description{

  This class extends
  \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}. It
  is a data frame and associated metadata (describing the columns of the
  data frame). The main purpose of this class is to contain alignment
  data in addition to the central information of
  \code{\linkS4class{AlignedRead}}.

}

\section{Objects from the Class}{

  Objects from the class are created by calls to the
  \code{\link{AlignedDataFrame}} function.

}
\section{Slots}{
  \describe{

    \item{\code{data}:}{Object of class \code{"data.frame"} containing
      the data. See \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} for details.}

    \item{\code{varMetadata}:}{Object of class \code{"data.frame"}
      describing columns of \code{data}.  See
      \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} for details.}

    \item{\code{dimLabels}:}{Object of class \code{character} describing
      the dimensions of the AnnotatedDataFrame. Used internally; see
      \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} for details.}

    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      describing the version of this object. Used internally; see
      \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} for details.}

  }
}
\section{Extends}{
Class \code{"\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}"}, directly.
Class \code{"\link[Biobase:class.Versioned]{Versioned}"}, by class "AnnotatedDataFrame", distance 2.
}
\section{Methods}{

  This class inherits methods \code{pData} (to retrieve the underlying
  data frame) and \code{varMetadata} (to retrieve the metadata) from
  \code{AnnotatedDataFrame}.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
}

\keyword{classes}
