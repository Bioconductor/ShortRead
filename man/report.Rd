\name{report}
\alias{report}

\title{Summarize quality assessment results into a report}

\description{

  This generic function summarizes results from evaluation of
  \code{\link{qa}} into a report. Available report formats vary
  depending on the data analysed.

}

\usage{
report(x, ..., dest = tempfile(), type = "pdf")
}

\arguments{
  \item{x}{An object returned by \code{\link{qa}}, usually derived from
    class \code{\linkS4class{.QA}}}

  \item{\dots}{Additional arguments used by specific methods. See
    specific methods for details.}

  \item{dest}{The output destination for the final report; destination
    files are overwritten if they already exist.}

  \item{type}{A text string defining the type of report; available
    report types depend on the type of object \code{x}.}
}

\value{
  This function is invoked for its side effect;  the return value is the
  name of the directory or file where the report was created.
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\linkS4class{SolexaExportQA}}
}

\examples{
showMethods("report")
}

\keyword{manip}
