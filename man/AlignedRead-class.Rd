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
\alias{coerce,AlignedRead,RangesList-method}
\alias{coerce,AlignedRead,RangedData-method}
\alias{coerce,AlignedRead,GRanges-method}
\alias{coerce,AlignedRead,GappedAlignments-method}
\alias{chromosome,AlignedRead-method}
\alias{position,AlignedRead-method}
\alias{strand,AlignedRead-method}
\alias{coverage,AlignedRead-method}
\alias{srrank,AlignedRead-method}
\alias{srorder,AlignedRead-method}
\alias{srduplicated,AlignedRead-method}
\alias{\%in\%,AlignedRead,RangesList-method}
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

    \item{\code{chromosome}}{Object of class \code{"factor"} the
      particular sequence within a set of target sequences
      (e.g. chromosomes in a genome assembly) to which each short read
      aligns.}

    \item{\code{position}}{Object of class \code{"integer"} the
      (base-pair) position in the genome to which the read is
      aligned. AlignedRead objects created by readAligned use 1-based
      indexing, with alignemnts reported in \sQuote{left-most}
      coordinates, as described in the vignette.}

    \item{\code{strand}}{Object of class \code{"factor"} the strand of
      the alignment.}

    \item{\code{alignQuality}}{Object of class \code{"numeric"}
      representing an alignment quality score.}

    \item{\code{alignData}}{Object of class \code{"AlignedDataFrame"}
      additional alignment information.}

    \item{\code{quality}}{Object of class \code{"BStringSet"}
      representing base-call read quality scores.}

    \item{\code{sread}}{Object of class \code{"DNAStringSet"} DNA
      sequence of the read.}

    \item{\code{id}}{Object of class \code{"BStringSet"} read
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

    \item{coerce}{
	  \code{signature(from = "PairwiseAlignedXStringSet", to = "AlignedRead")}:

	  \code{signature(from = "AlignedRead", to = "RangesList")}:
          \code{signature(from = "AlignedRead", to = "RangedData")}:
          \code{signature(from = "AlignedRead", to = "GRanges")}:
          \code{signature(from = "AlignedRead", to = "GappedAlignments")}:
      
      Invoke these methods with, e.g., \code{as(from, "AlignedRead")}
      to coerce objects of class \code{from} to class
      \code{"AlignedRead"}.
	  
      Coercion from AlignedRead to \link[IRanges]{RangesList},
          \link[IRanges]{RangedData} or
	  \link[GenomicRanges]{GRanges} assumes that
	  \code{position(from)} uses a \sQuote{leftmost} (see
	  \code{coverage} on this page) coordinate
	  system. Since \link[IRanges]{Ranges} objects cannot
	  store \code{NA} values, reads with \code{NA} in the
	  \code{position}, \code{width}, \code{chromosome} or (in the
	  case of GRanges) \code{strand} vectors are dropped.
  
	}

    \item{chromosome}{\code{signature(object = "AlignedRead")}: access the
      chromosome slot of \code{object}.}

    \item{position}{\code{signature(object = "AlignedRead")}: access the
      position slot of \code{object}.}

    \item{strand}{\code{signature(object = "AlignedRead")}: access the
      strand slot of \code{object}.}

    \item{coverage}{
      \code{signature(x = "AlignedRead", shift = 0L, width = NULL, weight = 1L, ..., coords = c("leftmost", "fiveprime"), extend=0L)}:

      Calculate coverage across reads present in \code{x}.

      \code{shift} must be either \code{0L} or a named integer vector
      with names including all \code{levels(chromosome(x))}.
      It specifies how the reads in \code{x} should be (horizontally)
      shifted \emph{before} the coverage is computed.

      \code{width} must be either \code{NULL} or a named vector of
      non-negative integers with names including all
      \code{levels(chromosome(x))}.  In the latter case, it specifies
      for each chromosome the end of the chromosome region over which
      coverage is to be calculated \emph{after} the reads have been
      shifted. Note that this region always starts at chromosome
      position 1.  If \code{width} is \code{NULL}, it ends at the
      rightmost chromosome position covered by at least one read.

      \code{weight} must be \code{1L} for now (weighting the reads is
      not supported yet, sorry).

      \code{coords} specifies the coordinate system used to record
      position. Both systems number base pairs from left to right on the
      5' strand. \code{leftmost} indicates the eland convention, where
      \code{position(x)} is the left-most (minimum) base pair,
      regardless of strand. \code{fiveprime} is the MAQ convention,
      where \code{position(x)} is the coordinate of the 5' end of the
      aligned read.

      \code{extend} indicates the number of base pairs to extend the
      read. Extension is in the 3' direction, measured from the
      3' end of the aligned read.

      The return value of \code{coverage} is a \code{SimpleRleList} object.
    }
    
	\item{\%in\%}{
	  \code{signature(x = "AlignedRead", table = "RangesList")}:

	  Return a length(x) logical vector indicating whether the
	  chromosome, position, and width of \code{x} overlap (see IRanges
	  \code{\link[IRanges:IntervalTree-class]{overlap}}) with ranges in
	  \code{table}. Reads for which \code{chromosome()},
	  \code{position()}, or \code{width()} return \code{NA} \emph{never}
	  overlap with \code{table}. This function assumes that positions	
	  are in \sQuote{leftmost} coordinates, as defined in
	  \code{coverage}.

	}
	
    \item{srorder}{\code{signature(x = "AlignedRead", ..., withSread=TRUE)}:}
    \item{srrank}{\code{signature(x = "AlignedRead", ..., withSread=TRUE)}:}
    \item{srsort}{\code{signature(x = "AlignedRead", ..., withSread=TRUE)}:}
    \item{srduplicated}{\code{signature(x = "AlignedRead", ..., withSread=TRUE)}:

      Order, rank, sort, and find duplicates in \code{AlignedRead}
      objects. Reads are sorted by \code{chromosome}, \code{strand},
      \code{position}, and then (if \code{withSread=TRUE}) \code{sread};
      less fine-grained sorting can be accomplished with, e.g.,
      \code{x[srorder(sread(x))]}. \code{srduplicated} behaves like
      \code{duplicated}, i.e., the first copy of a duplicate is
      \code{FALSE} while the remaining copies are \code{TRUE}.}

    \item{show}{\code{signature(object = "AlignedRead")}: provide a
      compact display of the \code{AlignedRead}  content.}

    \item{detail}{\code{signature(x = "AlignedRead")}: display
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
(aln <- readAligned(dirPath, 'out.aln.1.txt', type="MAQMapview"))
coverage(aln)[[1]]
cvg <- coverage(aln, shift=c(ChrA=10L))
## remove 0 coverage on left ends
ltrim0 <- function(x) {
    i <- !cumprod(runValue(x) == 0)
    Rle(runValue(x)[i], runLength(x)[i])
}
endoapply(cvg, ltrim0)
## demonstration of show() and detail() methods
show(aln)
detail(aln)
}

\keyword{classes}
