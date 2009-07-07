\name{readIntensities}

\alias{readIntensities}
\alias{readIntensities,character-method}

\title{Read Solexa 'int' and 'nse' files}

\description{

  \code{readIntensities} reads image `intensity' files (such as Solexa's
  \code{_int.txt} and (optionally) \code{_nse.txt}) in a directory into
  a single object.

}
\usage{
readIntensities(dirPath, pattern=character(0), ...)
}

\arguments{

  \item{dirPath}{Directory path or other object (e.g.,
    \code{\linkS4class{SolexaPath}} for which methods are defined.}

  \item{pattern}{A length 1 character vector representing a regular
    expression to be combined (using, e.g., \code{paste(pattern,
    intExtension, sep="")} with \code{intExtension} or
    \code{nseExtension} to match files to be summarized.}

  \item{\dots}{Additional arguments used by methods.}

}

\details{

  Additional methods are defined on specific classes, see, e.g.,
  \code{\linkS4class{SolexaPath}}.

  The \code{readIntensities,character-method} contains an argument
  \code{type} that determines how intensities are parsed. Use the
  \code{type} argument to \code{readIntensities,character-method}, as
  follows:

  \itemize{

    \item{\code{type="IparIntensity"}}{

      Intensities are read from Solexa \code{_pos.txt},
      \code{_int.txt.p}, \code{_nse.txt.p}-style file triplets. The
      signature for this method is

	  \code{dirPath, pattern=character(0), ...,
		type="IparIntensity",
        intExtension="_int.txt.p.gz",
        nseExtension="_nse.txt.p.gz",
        posExtension="_pos.txt",
        withVariability=TRUE, verbose=FALSE}

    }

    \item{\code{type="SolexaIntensity"}}{

      Intensities are read from Solexa \code{_int.txt} and
      \code{_nse.txt}-style files; see \code{\linkS4class{SolexaPath}}
      for details.  The signature for this method is

	  \code{dirPath, pattern=character(0), ...,
		type="SolexaIntensity",
        intExtension="_int.txt",
        nseExtension="_nse.txt",
        withVariability=TRUE, verbose=FALSE}

    }
  }

  Arguments to these methods are as follows:

  \itemize{

	\item{intExtension, nseExtension, posExtension}{\code{character(1)}
	  values pasted (with \code{sep=""}) to \code{pattern} to
	  identify different file sources.}

	\item{withVariability}{Include estimates of variability (i.e., from
	  parsing \code{_nse} files).}

	\item{verbose}{Report on progress when starting to read each file.}

  }
}

\value{
  An object derived from class \code{\linkS4class{Intensity}}.
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
fl <- system.file("extdata", package="ShortRead")
sp <- SolexaPath(fl)
int <- readIntensities(sp)
int
intensity(int)[1,,] 		# one read
intensity(int)[[1:2,,]]         # two reads, as 'array'
head(rowMeans(intensity(int)))	# treated as 'array'
head(pData(readInfo(int)))
}

\keyword{manip}
