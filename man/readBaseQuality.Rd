\name{readBaseQuality}

\alias{readBaseQuality}
\alias{readBaseQuality,character-method}

\title{Read short reads and their quality scores into R representations}

\description{

  \code{readBaseQuality} reads all base call files in a directory
  \code{dirPath} whose file name matches \code{seqPattern} and all quality
  score files whose name matches \code{prbPattern}, returning a compact
  internal representation of the sequences, and quality scores in the files.
  Methods read all files into a single R object.

}
\usage{

readBaseQuality(dirPath, ...)
\S4method{readBaseQuality}{character}(dirPath, seqPattern=character(0),
prbPattern=character(0), type=c("Solexa"), ...)

}

\arguments{

  \item{dirPath}{A character vector (or other object; see methods
    defined on this generic) giving the directory path (relative or
    absolute) of files to be input.}

  \item{seqPattern}{The (\code{\link{grep}}-style) pattern describing
    base call file names to be read. The default (\code{character(0)})
    results in (attempted) input of all files in the directory.}

  \item{prbPattern}{The (\code{\link{grep}}-style) pattern describing
    quality score file names to be read. The default (\code{character(0)})
    results in (attempted) input of all files in the directory.}

  \item{type}{The type of file to be parsed. Supported types include:
	\code{Solexa}: parse reads and their qualities from \code{_seq.txt}
	and \code{_prb.txt}-formatted files, respectively.}

  \item{...}{Additional arguments, perhaps used by methods.}

}

\value{

  A single R object (e.g., \code{\linkS4class{ShortReadQ}}) containing
  sequences and qualities of all files in \code{dirPath} matching
  \code{seqPattern} and \code{prbPattern} respectively. There is no
  guarantee of order in which files are read.

}

\seealso{

  A \code{\linkS4class{ShortReadQ}} object.
  
  \code{\link{readXStringColumns}}, \code{\link{readPrb}}

}

\author{
  Patrick Aboyoun <paboyoun@fhcrc.org>
}

\examples{
sp <- SolexaPath(system.file("extdata", package="ShortRead"))
readBaseQuality(sp, seqPattern="s_1.*_seq.txt", prbPattern="s_1.*_prb.txt")
}
\keyword{manip}
