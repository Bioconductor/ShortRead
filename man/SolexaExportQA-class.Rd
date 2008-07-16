\name{SolexaExportQA-class}
\docType{class}
\alias{SolexaExportQA}
\alias{SolexaExportQA-class}
\alias{report,SolexaExportQA-method}
\alias{show,SolexaExportQA-method}

\title{Quality assessment summaries from Solexa export files}

\description{

  This class contains a list-like structure with summary descriptions
  derived from visiting one or more Solexa \sQuote{export}
  files. 

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
  Accessor methods are inheritted from the \code{\linkS4class{SRList}}
  class.

  A additional methods defined on this class are:
  \describe{
    \item{report}{\code{signature(x="SolexaExportQA", ...,
        dest=paste(tempfile(), "pdf", sep="."), type="pdf")}: produces a
        pdf file summarizing the QA results.}

    \item{show}{\code{signature(object = "SolexaExportQA")}: Display an
      overview of the object contents.}
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\link{qa}}.
}
\examples{
showClass("SolexaExportQA")
}
\keyword{classes}
