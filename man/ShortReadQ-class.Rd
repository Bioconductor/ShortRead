\name{ShortReadQ-class}
\docType{class}

\alias{ShortReadQ-class}
\alias{ShortReadQ}
\alias{ShortReadQ,DNAStringSet,QualityScore,BStringSet-method}
\alias{ShortReadQ,DNAStringSet,QualityScore,missing-method}
\alias{ShortReadQ,DNAStringSet,BStringSet,BStringSet-method}
\alias{ShortReadQ,DNAStringSet,BStringSet,missing-method}
\alias{ShortReadQ,missing,missing,missing-method}
\alias{coerce,ShortReadQ,QualityScaledDNAStringSet-method}
\alias{writeFastq,ShortReadQ,character-method}
\alias{writeFastq,ShortReadQ,FastqFile-method}
\alias{alphabetByCycle,ShortReadQ-method}
\alias{alphabetScore,ShortReadQ-method}
\alias{[,ShortReadQ,ANY,ANY-method}
\alias{[,ShortReadQ,ANY,missing-method}
\alias{[,ShortReadQ,missing,ANY-method}
\alias{[,ShortReadQ,missing,missing-method}
\alias{[,ShortReadQ,missing,missing,ANY-method}
\alias{[,ShortReadQ,missing,ANY,ANY-method}
\alias{[,ShortReadQ,ANY,ANY,ANY-method}
\alias{[,ShortReadQ,ANY,missing,ANY-method}
\alias{append,ShortReadQ,ShortReadQ,missing-method}
\alias{narrow,ShortReadQ-method}
\alias{trimTailw,ShortReadQ-method}
\alias{trimTails,ShortReadQ-method}
\alias{qa,ShortReadQ-method}
\alias{detail,ShortReadQ-method}

\title{"ShortReadQ" class for short reads and their quality scores}

\description{

  This class provides a way to store and manipulate, in a coordinated
  fashion, the reads, identifiers, and quality scores of uniform-length
  short reads.

}

\section{Objects from the Class}{

  Objects from this class are the result of \code{\link{readFastq}}, or
  can be constructed from \code{DNAStringSet}, \code{QualityScore}, and
  \code{BStringSet} objects, as described below.

}
\section{Slots}{

  Slots \code{sread} and \code{id} are inherited from
  \code{\linkS4class{ShortRead}}. An additional slot defined in this
  class is:
  \describe{
    \item{\code{quality}:}{Object of class \code{"BStringSet"}
      representing a quality score (see \code{\link{readFastq}} for some
      discussion of quality score).}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ShortRead}"}, directly.
Class \code{"\linkS4class{.ShortReadBase}"}, by class "ShortRead", distance 2.
}
\section{Methods}{

  Constructors include:

  \describe{

    \item{ShortReadQ}{\code{signature(sread = "DNAStringSet", quality = "QualityScore", id = "BStringSet")}:}
    \item{ShortReadQ}{\code{signature(sread = "DNAStringSet", quality = "BStringSet", id = "BStringSet")}:

      Create a \code{ShortReadQ} object from reads, their quality
      scores, and identifiers. When \code{quality} is of class
      \code{BStringSet}, the type of encoded quality score is inferred
      from the letters used in the scores. The length of \code{id} and
      \code{quality} must match that of \code{sread}.}

    \item{ShortReadQ}{\code{signature(sread = "DNAStringSet", quality = "QualityScore", id = "missing")}:}
    \item{ShortReadQ}{\code{signature(sread = "DNAStringSet", quality = "BStringSet", id = "missing")}:

      Create a \code{ShortReadQ} object from reads and their quality
      scores, creating empty identifiers. When \code{quality} is of
      class \code{BStringSet}, the type of encoded quality score is
      inferred from the letters used in the scores. }

    \item{ShortReadQ}{\code{signature(sread = "missing", quality = "missing", id = "missing", ...)}:
      Create an empty \code{ShortReadQ} object.}

  }

  See \code{\link{accessors}} for additional functions to access slot
  content, and \code{\linkS4class{ShortRead}} for inherited
  methods. Additional methods include:

  \describe{

    \item{quality}{inherited from \code{signature(object = "ANY")}:
      access the quality slot of \code{object}.}


    \item{coerce}{\code{signature(from = "SFastqQuality", to = "QualityScaledDNAStringSet")}:
      
      (Use \code{as(from, "QualityScaledDNAStringSet")}) coerce objects
      of class \code{from} to class \code{to}, using the quality
      encoding implied by \code{quality(from)}. See
      \code{\linkS4class{QualityScore}} for supported quality
      classes and their coerced counterparts.}

    \item{writeFastq}{\code{signature(object = "ShortReadQ", file =
	"character", ...)}:}
    \item{writeFastq}{\code{signature(object = "ShortReadQ", file =
	"FastqFile", ...)}:
      Write \code{object} to \code{file} in fastq format. See
	\code{?\link{writeFastq}} for additional arguments \code{mode}
	and \code{full}.}

    \item{[}{\code{signature(x = "ShortReadQ", i = "ANY", j = "missing")}: 
      This method creates a new \code{ShortReadQ} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,ShortRead} do not provide additional functionality, but
      are present to limit inappropriate use.}

    \item{append}{\code{signature(x = "ShortReadQ", values = "ShortRead", length = "missing")}:
      append the \code{sread}, \code{quality} and \code{id} slots of
      \code{values} after the corresponding fields of \code{x}.}

    \item{narrow}{\code{signature(x = "ShortReadQ", start = NA, end =
      NA, width = NA, use.names = TRUE)}: narrow \code{sread} and
      \code{quality} so that sequences are between \code{start} and
      \code{end} bases, according to
      \code{\link[IRanges:intra-range-methods]{narrow}} in the \code{IRanges}
      package.}

    \item{trimTailw}{\code{signature(object="ShortReadQ", k="integer",
	a="character", halfwidth="integer", ..., ranges=FALSE)}:
      trim trailing nucleotides when a window of width 2 * halfwidth + 1
      contains \code{k} or more quality scores falling at or below
      \code{a}.}

    \item{trimTails}{\code{signature(object="ShortReadQ", k="integer",
	a="character", successive=FALSE, ..., ranges=FALSE)}: trim
      trailing nucleotides if \code{k} nucleotides fall below the
      quality encoded by \code{a}. If \code{successive=FALSE}, the k'th
      failing nucleotide and all subseqent nucleotides are trimmed. If
      \code{successive=TRUE}, failing nucleotides must occur
      successively; the sequence is trimmed from the first of the
      successive failing nucleotides.}

    \item{alphabetByCycle}{\code{signature(stringSet = "ShortReadQ")}:
      Apply \code{\link{alphabetByCycle}} to the \code{sread} component,
      the \code{quality} component, and the combination of these two
      components of \code{stringSet}, returning a list of matrices with
      three elements: \code{"sread"}, \code{"quality"}, and
      \code{"both"}.}

    \item{alphabetScore}{\code{signature(object = "ShortReadQ")}:
      See \code{\link{alphabetScore}} for details.}

    \item{qa}{\code{signature(dirPath = "ShortReadQ", lane="character",
	..., verbose=FALSE)}:
      Perform quality assessment on the \code{ShortReadQ} object using
      \code{lane} to identify the object and returning an instance of
      \code{\linkS4class{ShortReadQQA}}. See \code{\link{qa}}}

    \item{detail}{\code{signature(x = "ShortReadQ")}: display the
      first and last entries of each of \code{sread}, \code{id}, and
      \code{quality} entries of \code{object}.}
  }
}

\author{Martin Morgan}

\seealso{

  \code{\link{readFastq}} for creation of objects of this class from
  fastq-format files.

}
\examples{
showClass("ShortReadQ")
showMethods(class="ShortReadQ", where=getNamespace("ShortRead"),
            inherit=FALSE)
showMethods(class="ShortRead", where=getNamespace("ShortRead"),
            inherit=FALSE)

sp <- SolexaPath(system.file('extdata', package='ShortRead'))
rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
quality(rfq)
quality(trimTails(rfq, 2, "H", successive=TRUE))
}
\keyword{classes}
