\name{.QA-class}
\docType{class}
\alias{.QA-class}

\title{Virtual class for representing quality assessment results}
\description{

  Classes derived from \code{.QA-class} represent results of quality
  assurance analyses. Details of derived class structure are found on
  the help pages of the derived classes.

}

\section{Objects from the Class}{
  Objects from the class are created by ShortRead functions, in
  particular \code{\link{qa}}.
}

\section{Extends}{
  Class \code{"\linkS4class{.ShortReadBase}"}, directly.
}

\section{Methods}{

  There are no methods defined directly on the \code{QA} class; see
  derived class help pages for additional methods.

}

\seealso{
  \code{\linkS4class{SolexaExportQA}}.
}

\author{Martin Morgan <mtmmorgan@fhcrc.org>}

\examples{
getClass(".QA", where=getNamespace("ShortRead"))
}

\keyword{classes}
