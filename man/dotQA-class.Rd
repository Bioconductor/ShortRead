\name{.QA-class}
\docType{class}
\alias{.QA-class}
\alias{rbind,.QA-method}
\alias{show,.QA-method}

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
  
  Methods defined on this class include:

  \describe{

    \item{rbind}{\code{signature(...="list")}: rbind data frame objects
      in \code{...}. All objects of \code{...} must be of the same
      class; the return value is an instance of that class.}

    \item{show}{\code{signature(object = "SolexaExportQA")}: Display an
      overview of the object contents.}

  }
}

\seealso{
  Specific classes derived from \code{.QA}
}

\author{Martin Morgan <mtmmorgan@fhcrc.org>}

\examples{
getClass(".QA", where=getNamespace("ShortRead"))
}

\keyword{classes}
