\name{AlignedRead-class}
\docType{class}
\alias{AlignedRead-class}
\alias{[,AlignedRead,missing,missing-method}
\alias{[,AlignedRead,missing,ANY-method}
\alias{[,AlignedRead,ANY,ANY-method}
\alias{[,AlignedRead,ANY,missing-method}
\alias{detail,AlignedRead-method}
\alias{show,AlignedRead-method}

\title{"AlignedRead" class for aligned short reads}

\description{

  This class represents and manipulates reads and their genomic
  alignments. Alignment information includes genomic postion, strand,
  quality, and other data.

}

\section{Objects from the Class}{

  Objects of this class can be created rom a call to the
  \code{\link{AlignedRead}} constructor, or more typically by parsing
  appropriate files (e.g., \code{\link{readAligned}}).

}
\section{Slots}{
  \describe{

    \item{\code{chromosome}:}{Object of class \code{"factor"} containing
      the genome to which the read is aligned.}

    \item{\code{position}:}{Object of class \code{"integer"} the
    (base-pair) position in the genome to which the read is aligned.}

    \item{\code{strand}:}{Object of class \code{"factor"} the strand of
      the alignment.}

    \item{\code{alignQuality}:}{Object of class \code{"numeric"}
      represented an alignment quality score.}

    \item{\code{alignData}:}{Object of class \code{"AlignedDataFrame"}
      additional alignment information.}

    \item{\code{quality}:}{Object of class \code{"BStringSet"} base-call
      read quality scores.}

    \item{\code{sread}:}{Object of class \code{"DNAStringSet"} DNA
      sequence of the read.}

    \item{\code{id}:}{Object of class \code{"BStringSet"} read
      identifier.}

  }
}
\section{Extends}{
Class \code{"\linkS4class{ShortReadQ}"}, directly.
Class \code{"\linkS4class{ShortRead}"}, by class "ShortReadQ", distance 2.
Class \code{"\linkS4class{.ShortReadBase}"}, by class "ShortReadQ", distance 3.
}
\section{Methods}{
  \describe{

    \item{[}{\code{signature(x = "ShortReadQ", i = "ANY", j = "missing")}: 
      This method creates a new \code{ShortReadQ} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,ShortRead} do not provide additional functionality, but
      are present to limit inappropriate use.}

    \item{show}{\code{signature(object = "AlignedRead")}: provide a
      compact display of the \code{AlignedRead}  content.}

    \item{detail}{\code{signature(object = "AlignedRead")}: display
      \code{alignData} in more detail.}
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\link{readAligned}}
}
\examples{
showMethods(class="AlignedRead")
dirPath <- system.file('extdata', 'maq', package='ShortRead')
readAligned(dirPath, 'out.aln.1.txt', type="MAQMapview")
}
\keyword{classes}
