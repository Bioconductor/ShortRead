\name{ShortRead-class}
\docType{class}
\alias{ShortRead-class}
\alias{alphabetByCycle,ShortRead-method}
\alias{clean,ShortRead-method}
\alias{[,ShortRead,missing,missing-method}
\alias{[,ShortRead,missing,ANY-method}
\alias{[,ShortRead,ANY,ANY-method}
\alias{[,ShortRead,ANY,missing-method}
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

  Objects from this class are not normally created directly by the user.

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
  \describe{
    \item{[}{\code{signature(x = "ShortRead", i = "ANY", j = "missing")}: 
      This method creates a new \code{ShortRead} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,ShortRead} do not provide additional functionality, but
      are present to limit inappropriate use.}

    \item{length}{\code{signature(x = "ShortRead")}: returns a
      \code{integer(1)} vector describing the number of reads in this
      object.}

    \item{width}{\code{signature(x = "ShortRead")}: returns a
      \code{integer(1)} vector describing the number of characters in
      each read in this object.}

    \item{alphabetByCycle}{\code{signature(stringSet = "ShortRead")}:
      Apply \code{\link{alphabetByCycle}} to the \code{sread} component
      of \code{stringSet}, returning a matrix as described in
      \code{\link{alphabetByCycle}}.}

    \item{clean}{\code{signature(object="ShortRead")}: Remove all reads
      containing non-nucleotide (\code{"N", "-"}) symbols.}

    \item{show}{\code{signature(object = "ShortRead")}: provides a brief
      summary of the object, including its class, length and width.}

    \item{detail}{\code{signature(object = "ShortRead")}: provides a
      more extensive summary of this object, dispalying the first and
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
