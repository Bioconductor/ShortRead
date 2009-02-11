\name{FastqQA-class}
\docType{class}
\alias{FastqQA}
\alias{FastqQA-class}

\title{Quality assessment summaries from MAQ map files}

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
%   \describe{
%     \item{report}{\code{signature(x="FastqQA", ...,
%         dest=paste(tempfile(), "pdf", sep="."), type="pdf")}: produces a
%         pdf file summarizing the QA results.}
%   }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\link{qa}}.
}
\examples{
showClass("FastqQA")
}
\keyword{classes}
