\name{report}
\alias{report}
\alias{report,ANY-method}

\title{Summarize quality assessment results into a report}

\description{

  This generic function summarizes results from evaluation of
  \code{\link{qa}} into a report. Available report formats vary
  depending on the data analysed.

}

\usage{
report(x, ..., dest=tempfile(), type="html")
}

\arguments{
  \item{x}{An object returned by \code{\link{qa}}, usually derived from
    class \code{\linkS4class{.QA}}}

  \item{\dots}{Additional arguments used by specific methods. See
    specific methods for details.}

  \item{dest}{The output destination for the final report. For
    \code{type="html"} this is a directory; for (deprecated)
    \code{type="pdf"} this is a file.}

  \item{type}{A text string defining the type of report; available
    report types depend on the type of object \code{x}; usually this is
    \dQuote{html}.}
}

\details{
  The following methods are defined:
  \itemize{

    \item{\code{x="BowtieQA", ..., dest=tempfile(), type="html"}}{
	  Produce an HTML-based report from an object of class
          \code{\linkS4class{BowtieQA}}.}

    \item{\code{x="FastqQA", ..., dest=tempfile(), type="html"}}{
	  Produce an HTML-based report from an object of class
      \code{\linkS4class{FastqQA}}.
	}

    \item{\code{x="MAQMapQA", ..., dest=tempfile(), type="html"}}{
	  Produce an HTML-based report from an object of class
      \code{\linkS4class{MAQMapQA}}.
	}

    \item{\code{x="SolexaExportQA", ..., dest=tempfile(), type="html"}}{
	  Produce an HTML-based report from an object of class
      \code{\linkS4class{SolexaExportQA}}.
	}

    \item{\code{x="SolexaExportQA", ..., dest=tempfile(), type="pdf"}}{
	  (Deprecated) Produce an PDF report from an object of class
      \code{\linkS4class{SolexaExportQA}}.
	}

    \item{\code{x="SolexaPath", ..., dest=tempfile(), type="html"}}{
	  Produce an HTML report by first visiting all \code{_export.txt}
      files in the \code{analysisPath} directory of \code{x} to create a
      \code{SolexaExportQA} instance.
	}

    \item{\code{x="SolexaPath", ..., dest=tempfile(), type="pdf"}}{
	  (Deprecated) Produce an PDF report by first visiting all
      \code{_export.txt} files in the \code{analysisPath} directory of
      \code{x} to create a \code{SolexaExportQA} instance.
	}

    \item{
      \code{x="ANY", ..., dest=tempfile(), type="ANY"}
    }{This method is used internally}
  }
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
