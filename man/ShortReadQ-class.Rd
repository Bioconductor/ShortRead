\name{ShortReadQ-class}
\docType{class}

\alias{ShortReadQ-class}
\alias{ShortReadQ,DNAStringSet,QualityScore,BStringSet-method}
\alias{ShortReadQ,DNAStringSet,QualityScore,missing-method}
\alias{ShortReadQ,ShortReadQ,missing,missing-method}
\alias{ShortReadQ,missing,missing,missing-method}
\alias{writeFastq,ShortReadQ-method}
\alias{alphabetByCycle,ShortReadQ-method}
\alias{alphabetScore,ShortReadQ-method}
\alias{[,ShortReadQ,missing,missing,ANY-method}
\alias{[,ShortReadQ,missing,ANY,ANY-method}
\alias{[,ShortReadQ,ANY,ANY,ANY-method}
\alias{[,ShortReadQ,ANY,missing,ANY-method}
\alias{append,ShortReadQ,ShortReadQ,missing-method}
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

    \item{ShortReadQ}{\code{signature(sread = "DNAStringSet", quality = "QualityScore", id = "BStringSet")}:
      Create a \code{ShortReadQ} object from reads, their quality
      scores, and identifiers. The length of \code{id} and
      \code{quality}must match that of \code{sread}.}

    \item{ShortReadQ}{\code{signature(sread = "DNAStringSet", quality = "QualityScore", id = "missing")}:
      Create a \code{ShortReadQ} object from reads and their quality
      scores, creating empty identifiers.}

    \item{ShortReadQ}{\code{signature(sread = "ShortReadQ", quality = "missing", id = "missing", ..., start=NA, end=NA)}:
      Create a \code{ShortReadQ} object from a previously existing
      \code{ShortReadQ} object, using the vectors \code{start},
      \code{end} (typically of length 1 or the length of \code{sread})
      to select the start and end of each read and quality
      score. \code{start} and \code{end} follow the semanitics of
      corresponding arguments to \code{DNAStringSet}}

    \item{ShortReadQ}{\code{signature(sread = "missing", quality = "missing", id = "missing", ...)}:
      Create an empty \code{ShortReadQ} object.}

  }

  See \code{\link{accessors}} for additional functions to access slot
  content, and \code{\linkS4class{ShortRead}} for inherited
  methods. Additional methods include:

  \describe{

    \item{writeFastq}{\code{signature(object = "ShortReadQ", file = "character", mode="character", ...)}:
      Write \code{object} to \code{file} in fastq format. \code{mode}
      defaults to \sQuote{w}. This creates a new file, or fails if
      \code{file} already exists. Use \code{mode="a"} to append to an
      existing file. \code{file} is expanded using \code{path.expand}.}

    \item{[}{\code{signature(x = "ShortReadQ", i = "ANY", j = "missing")}: 
      This method creates a new \code{ShortReadQ} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,ShortRead} do not provide additional functionality, but
      are present to limit inappropriate use.}

    \item{append}{\code{signature(x = "ShortRead", values = "ShortRead", length = "missing")}:
      append the \code{sread}, \code{quality} and \code{id} slots of
      \code{values} after the corresponding fields of \code{x}.}

    \item{alphabetByCycle}{\code{signature(stringSet = "ShortReadQ")}:
      Apply \code{\link{alphabetByCycle}} to the \code{sread} component,
      the \code{quality} component, and the combination of these two
      components of \code{stringSet}, returning a list of matrices with
      three elements: \code{"sread"}, \code{"quality"}, and
      \code{"both"}.}

    \item{alphabetScore}{\code{signature(object = "ShortReadQ")}:
      See \code{\link{alphabetScore}} for details.}

    \item{detail}{\code{signature(object = "ShortReadQ")}: display the
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
showMethods(class="ShortReadQ", inherit=FALSE)
showMethods(class="ShortRead", inherit=FALSE)
}
\keyword{classes}
