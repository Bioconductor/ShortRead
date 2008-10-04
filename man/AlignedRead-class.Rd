\name{AlignedRead-class}
\docType{class}
\alias{AlignedRead-class}
\alias{[,AlignedRead,missing,missing,ANY-method}
\alias{[,AlignedRead,missing,ANY,ANY-method}
\alias{[,AlignedRead,ANY,ANY,ANY-method}
\alias{[,AlignedRead,ANY,missing,ANY-method}
\alias{strand,AlignedRead-method}
\alias{detail,AlignedRead-method}
\alias{srrank,AlignedRead-method}
\alias{srorder,AlignedRead-method}
\alias{srduplicated,AlignedRead-method}
\alias{show,AlignedRead-method}

\title{"AlignedRead" class for aligned short reads}

\description{

  This class represents and manipulates reads and their genomic
  alignments. Alignment information includes genomic position, strand,
  quality, and other data.

}

\section{Objects from the Class}{

  Objects of this class can be created from a call to the
  \code{\link{AlignedRead}} constructor, or more typically by parsing
  appropriate files (e.g., \code{\link{readAligned}}).

}
\section{Slots}{
  \describe{

    \item{\code{chromosome}:}{Object of class \code{"factor"} the
      particular sequence within a set of target sequences
      (e.g. chromosomes in a genome assembly) to which each short read
      aligns.}

    \item{\code{position}:}{Object of class \code{"integer"} the
    (base-pair) position in the genome to which the read is aligned.}

    \item{\code{strand}:}{Object of class \code{"factor"} the strand of
      the alignment.}

    \item{\code{alignQuality}:}{Object of class \code{"numeric"}
      representing an alignment quality score.}

    \item{\code{alignData}:}{Object of class \code{"AlignedDataFrame"}
      additional alignment information.}

    \item{\code{quality}:}{Object of class \code{"BStringSet"}
      representing base-call read quality scores.}

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

  See \code{\link{accessors}} for additional functions to access slot
  content, and \code{\linkS4class{ShortReadQ}},
  \code{\linkS4class{ShortRead}} for inherited methods. Additional
  methods include:

  \describe{

    \item{[}{\code{signature(x = "ShortReadQ", i = "ANY", j = "missing")}: 
      This method creates a new \code{ShortReadQ} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,ShortRead} do not provide additional functionality, but
      are present to limit inappropriate use.}

    \item{strand}{\code{signature(object = "AlignedRead")}: access the
      strand slot of \code{object}.}

    \item{srorder}{\code{\signature{x = "AlignedRead"}}:}
    \item{srrank}{\code{\signature{x = "AlignedRead"}}:}
    \item{srsort}{\code{\signature{x = "AlignedRead"}}:}
    \item{srduplicated}{\code{\signature{x = "AlignedRead"}}:
      Order, rank, sort, and find duplicates in \code{AlignedRead}
      objects. Reads are sorted by \code{chromosome}, \code{strand},
      \code{position}, and then \code{sread}; less fine-grained sorting
      can be accomplished with, e.g., \code{x[srorder(sread(x))]}.}

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
