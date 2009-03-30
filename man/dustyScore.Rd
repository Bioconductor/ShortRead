\name{dustyScore}

\alias{dustyScore}
\alias{dustyScore,DNAStringSet-method}
\alias{dustyScore,ShortRead-method}

\title{Summarize low-complexity sequences}

\description{

  \code{dustyScore} identifies low-complexity sequences, in a manner
  inspired by the \code{dust} implementation in \code{BLAST}.

}
\usage{

dustyScore(x, ...)

}

\arguments{

  \item{x}{A \code{DNAStringSet} object, or object derived from
    \code{ShortRead}, containing a collection of reads to be
    summarized.}

  \item{...}{Additional arguments, not currently used.}

}

\details{

  The following methods are defined:
  \describe{

    \item{dustyScore}{\code{signature(x = "DNAStringSet")}: operating on
      an object derived from class \code{DNAStringSet}.}

    \item{dustyScore}{\code{signature(x = "ShortRead")}: operating on
      the \code{sread} of an object derived from class
      \code{ShortRead}.}

  }

  The dust-like calculations used here are as implemented at
  \url{https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2009-February/000170.html}.
  Scores range from 0 (all triplets unique) to the square of the width
  of the longest sequence (poly-A, -C, -G, or -T). 

}

\value{

  A vector of numeric scores, with length equal to the length of
  \code{x}.

}

\references{

  Morgulis, Getz, Schaffer and Agarwala, 2006. WindowMasker: window-based
  masker for sequenced genomes, Bioinformatics 22: 134-141.

}

\seealso{

  The WindowMasker supplement defining \code{dust}
  \url{ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/windowmasker/windowmasker_suppl.pdf}

}

\author{Herve Pages (code); Martin Morgan}

\examples{
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
range(dustyScore(rfq))
}
\keyword{manip}
