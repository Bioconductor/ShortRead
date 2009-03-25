\name{ShortRead-class}
\docType{class}
\alias{ShortRead-class}
\alias{ShortRead}
\alias{ShortRead,DNAStringSet,BStringSet-method}
\alias{ShortRead,DNAStringSet,missing-method}
\alias{ShortRead,missing,missing-method}
\alias{srrank,ShortRead-method}
\alias{srorder,ShortRead-method}
\alias{srsort,ShortRead-method}
\alias{srduplicated,ShortRead-method}
\alias{srdistance,ShortRead,ANY-method}
\alias{trimLRPatterns,ShortRead-method}
\alias{alphabetByCycle,ShortRead-method}
\alias{tables,ShortRead-method}
\alias{clean,ShortRead-method}
\alias{[,ShortRead,ANY,ANY-method}
\alias{[,ShortRead,ANY,missing-method}
\alias{[,ShortRead,missing,ANY-method}
\alias{[,ShortRead,missing,missing-method}
\alias{[,ShortRead,missing,missing,ANY-method}
\alias{[,ShortRead,missing,ANY,ANY-method}
\alias{[,ShortRead,ANY,ANY,ANY-method}
\alias{[,ShortRead,ANY,missing,ANY-method}
\alias{append,ShortRead,ShortRead,missing-method}
\alias{narrow,ShortRead-method}
\alias{detail,ShortRead-method}
\alias{length,ShortRead-method}
\alias{show,ShortRead-method}
\alias{width,ShortRead-method}

\title{"ShortRead" class for short reads}

\description{

  This class provides a way to store and manipulate, in a coordinated
  fashion, uniform-length short reads and their identifiers.

}

\section{Objects from the Class}{

  Objects from this class are created by \code{readFasta}, or by calls
  to the constructor \code{ShortRead}, as outlined below.

}
\section{Slots}{
  \describe{
    \item{\code{sread}:}{Object of class \code{"DNAStringSet"}
      containing IUPAC-standard, uniform-length DNA strings represent
      short sequence reads.}
    \item{\code{id}:}{Object of class \code{"BStringSet"} containing
      identifiers, one for each short read.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{.ShortReadBase}"}, directly.
}
\section{Methods}{

  Constructors include:

  \describe{

    \item{ShortRead}{\code{signature(sread = "DNAStringSet", id = "BStringSet")}:
      Create a \code{ShortRead} object from reads and their
      identifiers. The length of \code{id} must match that of \code{sread}.}

    \item{ShortRead}{\code{signature(sread = "DNAStringSet", id = "missing")}:
      Create a \code{ShortRead} object from reads, creating empty identifiers.}

    \item{ShortRead}{\code{signature(sread = "missing", id = "missing", ...)}:
      Create an empty \code{ShortRead} object.}

  }

  See \code{\link{accessors}} for slot accessor functions.

  \describe{
    \item{[}{\code{signature(x = "ShortRead", i = "ANY", j = "missing")}: 
      This method creates a new \code{ShortRead} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,ShortRead} do not provide additional functionality, but
      are present to limit inappropriate use.}

    \item{append}{\code{signature(x = "ShortRead", values = "ShortRead", length = "missing")}:
      append the \code{sread} and \code{id} slots of \code{values} after
      the corresponding fields of \code{x}.}

    \item{narrow}{\code{signature(x = "ShortRead", start = NA, end = NA, width = NA, use.names = TRUE)}:
      \sQuote{narrow} \code{sread} so that sequences are between
      \code{start} and \code{end} bases, according to
      \code{\link[IRanges:IRanges-utils]{narrow}} in the \code{IRanges}
      package.
    }

    \item{length}{\code{signature(x = "ShortRead")}: returns a
      \code{integer(1)} vector describing the number of reads in this
      object.}

    \item{width}{\code{signature(x = "ShortRead")}: returns an
      \code{integer()} vector of the widths of each read in this
      object.}

    \item{srorder}{\code{signature(x = "ShortRead")}:}
    \item{srrank}{\code{signature(x = "ShortRead")}:}
    \item{srsort}{\code{signature(x = "ShortRead")}:}
    \item{srduplicated}{\code{signature(x = "ShortRead")}:
      Order, rank, sort, and find duplicates in \code{ShortRead} objects
      based on \code{sread(x)}, analogous to the corresponding functions
      \code{order}, \code{rank}, \code{sort}, and \code{duplicated},
      ordering nucleotides in the order \code{ACGT}.}

    \item{srdistance}{\code{signature(pattern="ShortRead", subject="ANY")}:
      Find the edit distance between each read in \code{pattern} and the
      (short) sequences in \code{subject}. See \code{\link{srdistance}}
      for allowable values for \code{subject}, and for additional
      details.}

    \item{trimLRPatterns}{\code{signature(Lpattern = "", Rpattern = "", subject = "ShortRead", max.Lmismatch = 0, max.Rmismatch = 0, with.Lindels = FALSE, with.Rindels = FALSE, Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)}:

      Remove left and / or right flanking patterns from
      \code{sread(subject)}, as described in
      \code{\link[Biostrings:trimLRPatterns]{trimLRPatterns}}. Classes
      derived from \code{ShortRead} (e.g., \code{\link{ShortReadQ}},
      \code{\link{AlignedRead}}) have corresponding base quality scores
      trimmed, too. A user-supplied argument \code{ranges} is ignored by
      this method; the class of the return object is the same as the
      class of \code{subject}.}

    \item{alphabetByCycle}{\code{signature(stringSet = "ShortRead")}:
      Apply \code{\link{alphabetByCycle}} to the \code{sread} component
      of \code{stringSet}, returning a matrix as described in
      \code{\link{alphabetByCycle}}.}

    \item{tables}{\code{signature(x= "ShortRead", n = 50)}:
      Apply \code{\link{tables}} to the \code{sread} component of
      \code{x}, returning a list summarizing frequency of reads in
      \code{x}.}

    \item{clean}{\code{signature(object="ShortRead")}: Remove all reads
      containing non-nucleotide (\code{"N", "-"}) symbols.}

    \item{show}{\code{signature(object = "ShortRead")}: provides a brief
      summary of the object, including its class, length and width.}

    \item{detail}{\code{signature(object = "ShortRead")}: provides a
      more extensive summary of this object, displaying the first and
      last entries of \code{sread} and \code{id}.}
  }
}
\author{Martin Morgan}
\seealso{
  \code{\linkS4class{ShortReadQ}}
}
\examples{
showClass("ShortRead")
showMethods(class="ShortRead")
}
\keyword{classes}
