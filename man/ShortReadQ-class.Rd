\name{ShortReadQ-class}
\docType{class}

\alias{ShortReadQ-class}
\alias{alphabetByCycle,ShortReadQ-method}
\alias{[,ShortReadQ,missing,missing-method}
\alias{[,ShortReadQ,missing,ANY-method}
\alias{[,ShortReadQ,ANY,ANY-method}
\alias{[,ShortReadQ,ANY,missing-method}
\alias{detail,ShortReadQ-method}

\title{"ShortReadQ" class for short reads and their quality scores}

\description{

  This class provides a way to store and manipulate, in a coordinated
  fashion, the reads, identifiers, and quality scores of uniform-length
  short reads.

}

\section{Objects from the Class}{

  Objects from this class are not normally created directly by the user,
  but are instead the result of \code{\link{readFastq}}.

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
  This class inherits methods from \code{\linkS4class{ShortRead}}, and
  has the following additional methods defined:
  \describe{

    \item{[}{\code{signature(x = "ShortReadQ", i = "ANY", j = "missing")}: 
      This method creates a new \code{ShortReadQ} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,ShortRead} do not provide additional functionality, but
      are present to limit inappropriate use.}

    \item{alphabetByCycle}{\code{signature(stringSet = "ShortReadQ")}:
      Apply \code{\link{alphabetByCycle}} to the \code{sread} component,
      the \code{quality} component, and the combination of these two
      components of \code{stringSet}, returning a list of matrices with
      three elements:  \code{"sread"}, \code{"quality"}, and \code{"both"}.}

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
