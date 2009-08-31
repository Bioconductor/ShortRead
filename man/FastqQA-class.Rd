\name{ShortReadQA-class}
\docType{class}
\alias{FastqQA}
\alias{FastqQA-class}
\alias{ShortReadQQA-class}
\alias{report,FastqQA-method}

\title{Quality assessment of fastq files and ShortReadQ objects}

\description{

  These classes contains a list-like structure with summary descriptions
  derived from visiting one or more fastq files, or from a
  \code{\linkS4class{ShortReadQ}} object.

}

\section{Objects from the Class}{
  
  Objects of the class are usually produced by a \code{\link{qa}}
  method.

}
\section{Slots}{
  \describe{
    \item{\code{.srlist}:}{Object of class \code{"list"}, containing
      data frames or lists of data frames summarizing the results of
      qa.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SRList}"}, directly.
Class \code{"\linkS4class{.QA}"}, directly.
Class \code{"\linkS4class{.SRUtil}"}, by class "SRList", distance 2.
Class \code{"\linkS4class{.ShortReadBase}"}, by class ".QA", distance 2.
}

\section{Methods}{
  Accessor methods are inherited from the \code{\linkS4class{SRList}}
  class.

  Additional methods defined on this class are:
  \describe{

    \item{report}{\code{signature(x="FastqQA", ..., dest=tempfile(), type="html")}:
	  produces HTML files summarizing QA results. \code{dest} should be
	  a directory.}

    \item{report}{\code{signature(x="ShortReadQA", ..., dest=tempfile(), type="html")}:
	  produces HTML files summarizing QA results. \code{dest} should be
	  a directory.}

  }

}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\link{qa}}.
}
\examples{
showClass("FastqQA")
}
\keyword{classes}
