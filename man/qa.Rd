\name{qa}
\alias{qa}
\alias{qa,character-method}
\alias{qa,list-method}

\title{Perform quality assessment on short reads}

\description{

  This function is a common interface to quality assessment functions
  available in \code{ShortRead}. Results from this function may be
  displayed in brief, or integrated into reports using, e.g.,
  \code{\link{report}}.

}
\usage{
qa(dirPath, ...)
\S4method{qa}{character}(dirPath, pattern=character(0), 
    type=c("SolexaExport", "SolexaRealign", "Bowtie", "MAQMap", 
           "MAQMapShort", "fastq"),
    ...)
\S4method{qa}{list}(dirPath, ...)
}

\arguments{

  \item{dirPath}{A character vector or other object (e.g.,
    \code{\link{SolexaPath}}; see \code{showMethods}, below) locating
    the data for which quality assessment is to be performed. See help
    pages for defined methods (by evaluating the example code, below)
    for details of available methods.}

  \item{pattern}{A character vector limiting the files in \code{dirPath}
	to be processed, as with \code{\link{list.files}}. Care should be
	taken to specify pattern to avoid reading unintended files.}

  \item{type}{The type of file being parsed; must be a character vector
	of length 1, selected from one of the types enumerated in the
	parameter.}

  \item{\dots}{Additional arguments used by methods.}
}

\details{
  The following methods are defined, in addition to those on S4 formal
  classes documented elsewhere:

  \describe{

    \item{\code{qa,character-method}}{

      Quality assessment is performed on all files in directory
      \code{dirPath} whose file name matches \code{pattern}. The type of
      analysis performed is based on the \code{type} argument. Use
      \code{SolexaExport} when all files matching \code{pattern} are
      Solexa \code{_export.txt} files. Use \code{SolexaRealign} for
      Solexa \code{_realign.txt} files. Use \code{Bowtie} for Bowtie
      files. Use \code{MAQMapShort} for MAQ \code{map} files produced by
      MAQ versions below 0.70 and \code{MAQMap} for more recent output.
      Use \code{fastq} for collections of fastq-format files. Quality
      assessment details vary depending on data source.

    }
	
    \item{\code{qa,list-method}}{

	  \code{dirPath} is a list of objects, all of the same class and
	  typically derived from \code{ShortReadQ}, on which quality
	  assessment is performed. All elements of the list must have names,
	  and these should be unique.

    }
  }
}

\value{
  An object derived from class \code{\linkS4class{.QA}}
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\linkS4class{.QA}},
  \code{\linkS4class{SolexaExportQA}}
  \code{\linkS4class{MAQMapQA}}
  \code{\linkS4class{FastqQA}}
}

\examples{
showMethods("qa", where=getNamespace("ShortRead"))
}

\keyword{manip}
