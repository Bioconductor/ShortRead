\name{qa}
\alias{qa}

\title{Perform quality assessment on short reads}

\description{

  This function is a common interface to quality assessment functions
  available in \code{ShortRead}. Results from this function may be
  displayed in brief, or integrated into reports using, e.g.,
  \code{\link{report}}.

}
\usage{
qa(dirPath, ...)
}

\arguments{

  \item{dirPath}{A character vector or other object locating the data
    for which quality assessment is to be performed. See help pages for
    defined methods (by evaluating the example code, below) for details
    of available methods.}

  \item{\dots}{Additional arguments used by methods.}
}

\details{
  The following general method is defined, in addition to those for
  specific objects:
  \itemize{
    \item{
      \code{dirPath="character", pattern="character", type="SolexaExport", ...}}{

      Quality assessment is performed on all files in directory
      \code{dirPath} whose file name matches \code{pattern}. The type of
      analysis performed is based on the \code{type} argument.
    }
  }
}

\value{
  An object derived from class \code{\linkS4class{.QA}}
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\linkS4class{.QA}}, \code{\linkS4class{SolexaExportQA}}
}

\examples{
showMethods("qa")
}

\keyword{manip}
