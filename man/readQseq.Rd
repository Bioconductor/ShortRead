\name{readQseq}

\alias{readQseq}
\alias{readQseq,character-method}

\title{Read Solexa qseq files as fastq-style quality scores}

\description{

  \code{readQseq} reads all files matching \code{pattern} in a directory
  into a single \code{\linkS4class{ShortReadQ}}-class object.
  Information on machine, lane, tile, x, and y coordinates, filtering
  status, and read number are not returned (although filtering status
  can be used to selectively include reads as descrbied below).

}
\usage{
readQseq(dirPath, pattern = character(0), ...,
         as=c("ShortReadQ", "XDataFrame"),
         filtered=FALSE,
         verbose=FALSE)
}

\arguments{

  \item{dirPath}{Directory path or other object (e.g.,
    \code{\linkS4class{SolexaPath}}) for which methods are defined.}

  \item{pattern}{Regular expression matching names of \code{_qseq} files
    to be summarized.}

  \item{\dots}{Additional argument, passed to I/O functions.}

  \item{as}{\code{character(1)} indicating the class of the return
	type.}

  \item{filtered}{\code{logical(1)} indicating whether to include only
	those reads passing Solexa filtering?}

  \item{verbose}{\code{logical(1)} indicating whether to report on
	progress during evaluation.}

}

\value{
  An object of class \code{\linkS4class{ShortReadQ}}.
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
fl <- system.file("extdata", package="ShortRead")
sp <- SolexaPath(fl)
readQseq(sp)
}

\keyword{manip}
