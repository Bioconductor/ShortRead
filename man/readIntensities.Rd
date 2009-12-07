\name{readIntensities}

\alias{readIntensities}
\alias{readIntensities,character-method}

\title{Read Illumina image intensity files}

\description{

  \code{readIntensities} reads image `intensity' files (such as
  Illumina's \code{_int.txt} and (optionally) \code{_nse.txt}) into a
  single object.

}
\usage{
readIntensities(dirPath, pattern=character(0), ...)
}

\arguments{

  \item{dirPath}{Directory path or other object (e.g.,
    \code{\linkS4class{SolexaPath}}) for which methods are defined.}

  \item{pattern}{A length 1 character vector representing a regular
    expression to be combined with \code{dirPath}, as described below,
    to match files to be summarized.}

  \item{\dots}{Additional arguments used by methods.}

}

\details{

  Additional methods are defined on specific classes, see, e.g.,
  \code{\linkS4class{SolexaPath}}.

  The \code{readIntensities,character-method} contains an argument
  \code{type} that determines how intensities are parsed. Use the
  \code{type} argument to \code{readIntensities,character-method}, as
  described below. All \code{readIntensities,character} methods accepts
  the folling arguments:

  \describe{
	
	\item{withVariability:}{Include estimates of variability (i.e., from
	  parsing \code{_nse} files).}

	\item{verbose:}{Report on progress when starting to read each file.}

  }

  The supported types and their signatures are:

  \describe{
	
    \item{\code{type="RtaIntensity"}}{

      Intensities are read from Illumina \code{_cif.txt} and
      \code{_cnf.txt}-style files.  The signature for this method is

	  \code{dirPath, pattern=character(0), ..., type="RtaIntensity",
        lane=integer(0), cycles=integer(0), cycleIteration=1L,
        tiles=integer(0),
        laneName=sprintf("L%.3d", lane),
        cycleNames=sprintf("C%d.%d", cycles, cycleIteration),
        tileNames=sprintf("s_%d_%d", lane, tiles),
        posNames=sprintf("s_%d_%.4d_pos.txt", lane, tiles),
        withVariability=TRUE, verbose=FALSE}

	  \describe{

		\item{lane:}{\code{integer(1)} identifying the lane in which
		  cycles and tiles are to be processed.}

		\item{cycles:}{\code{integer()} enumerating cycles to be
		  processed.}

		\item{cycleIteration:}{\code{integer(1)} identifying the
		  iteration of the base caller to be summarized}

		\item{tiles:}{\code{integer()} enumerating tile numbers to be
		  summarized.}

		\item{laneName, cycleNames, tileNames,
		  posNames:}{\code{character()} vectors identifying the lane and
		  cycle directories, and the \sQuote{pos} and tile file names
		  (excluding the \sQuote{.cif} or \sQuote{.cnf} extension) to be
		  processed.}

	  }

	  The \code{dirPath} and \code{pattern} arguments are combined as
	  \code{list.files(dirPath, pattern)}, and must identify a single
	  directory. Most uses of this function will focus on a single tile
	  (specified with, e.g., \code{tiles=1L}); the \code{laneName},
	  \code{cycleNames}, \code{tileNames}, and \code{posNames}
	  parameters are designed to work with the default Illumina pipeline
	  and do not normally need to be specified.

    }

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

	  Files to be parsed are determined as, e.g., \code{paste(pattern,
	  intExtension, sep="")}.

    }

    \item{\code{type="SolexaIntensity"}}{

      Intensities are read from Solexa \code{_int.txt} and
      \code{_nse.txt}-style files.  The signature for this method is

	  \code{dirPath, pattern=character(0), ...,
		type="SolexaIntensity",
        intExtension="_int.txt",
        nseExtension="_nse.txt",
        withVariability=TRUE, verbose=FALSE}

	  Files to be parsed are determined as, e.g., \code{paste(pattern,
	  intExtension, sep="")}.
    }

  }

}

\value{
  An object derived from class \code{\linkS4class{Intensity}}.
}

\author{
  Martin Morgan <mtmorgan@fhcrc.org>,
  Michael Muratet <mmuratet@hudsonalpha.org> (RTA).
}

\examples{
fl <- system.file("extdata", package="ShortRead")
sp <- SolexaPath(fl)
int <- readIntensities(sp)
int
intensity(int)[1,,] 		# one read
intensity(int)[[1:2,,]]         # two reads, as 'array'
head(rowMeans(intensity(int)))	# treated as 'array'
head(pData(readInfo(int)))

\dontrun{## RTA Lane 2, cycles 1:80, cycle iteration 1, tile 3
int <- readIntensities("Data/Intensities", type="RtaIntensity",
                       lane=2, cycles=1:80, tiles=3)
}
}

\keyword{manip}
