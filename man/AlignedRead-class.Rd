\name{AlignedRead-class}
\docType{class}
\alias{AlignedRead-class}
\alias{[,AlignedRead,ANY,ANY-method}
\alias{[,AlignedRead,ANY,missing-method}
\alias{[,AlignedRead,missing,ANY-method}
\alias{[,AlignedRead,missing,missing-method}
\alias{[,AlignedRead,missing,missing,ANY-method}
\alias{[,AlignedRead,missing,ANY,ANY-method}
\alias{[,AlignedRead,ANY,ANY,ANY-method}
\alias{[,AlignedRead,ANY,missing,ANY-method}
\alias{append,AlignedRead,AlignedRead,missing-method}
\alias{coerce,PairwiseAlignedXStringSet,AlignedRead-method}
\alias{strand,AlignedRead-method}
\alias{coverage,AlignedRead-method}
\alias{srrank,AlignedRead-method}
\alias{srorder,AlignedRead-method}
\alias{srduplicated,AlignedRead-method}
\alias{detail,AlignedRead-method}
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

    \item{[}{\code{signature(x = "AlignedRead", i = "ANY", j = "missing")}: 
      This method creates a new \code{AlignedRead} object containing only
      those reads indexed by \code{i}. \code{chromosome} is recoded to
      contain only those levels in the new subset.}

    \item{append}{\code{signature(x = "AlignedRead", values = "AlignedRead", length = "missing")}:
      append \code{values} after \code{x}. \code{chromosome} and
      \code{strand} must be factors with the same levels. See methods
      for \code{ShortReadQ}, \code{AlignedDataFrame} for details of how
      these components of \code{x} and \code{y} are appended.}

    \item{coerce}{\code{signature(from = "PairwiseAlignedXStringSet", to = "AlignedRead")}:
      
      (Invoke this method with, \code{as(from, "AlignedRead")}) coerce
      objects of class \code{from} to class \code{to}.}

    \item{strand}{\code{signature(object = "AlignedRead")}: access the
      strand slot of \code{object}.}

    \item{coverage}{\code{signature(x = "AlignedRead", start = NA, end = NA, ..., coords=c("leftmost", "fiveprime"), extend=0L)}:

      Calculate coverage across reads present in \code{x}. 

      \code{start} and \code{end} are regions (e.g., of chromosomes)
      over which coverage is to be calculated. If provided, these are
      length 1 integers or \emph{named} integer vectors of length
      greater than 1. If named integer vectors, the names must match
      \code{levels(chromosome(x))}. If omitted, coverage is calculated
      over the range of values spanned by the reads in \code{x}

      \code{coords} specifies the coordinate system used to record
      position. Both systems number base pairs from left to right on the
      5' strand. \code{leftmost} indicates the eland convention,
      where \code{position(x)} is the left-most (minimum) base pair,
      regardless of strand. \code{fiveprime} is the MAQ convetion, where
      \code{position(x)} is the coordinate of the 5' end of the
      aligned read.

      \code{extend} indicates the number of base pairs to extend the
      read. Extension is in the 3' direction, measured from the
      3' end of the aligned read.

    }
    

    \item{srorder}{\code{signature{x = "AlignedRead"}}:}
    \item{srrank}{\code{signature{x = "AlignedRead"}}:}
    \item{srsort}{\code{signature{x = "AlignedRead"}}:}
    \item{srduplicated}{\code{signature{x = "AlignedRead"}}:
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
showMethods(class="AlignedRead", where=getNamespace("ShortRead"))
dirPath <- system.file('extdata', 'maq', package='ShortRead')
readAligned(dirPath, 'out.aln.1.txt', type="MAQMapview")
}
\keyword{classes}
