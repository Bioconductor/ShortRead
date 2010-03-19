\name{QualityScore-class}
\docType{class}
\alias{QualityScore-class}
\alias{NumericQuality-class}
\alias{IntegerQuality-class}
\alias{MatrixQuality-class}
\alias{FastqQuality-class}
\alias{SFastqQuality-class}
% QualityScore
\alias{[,QualityScore,ANY,missing-method}
\alias{[,QualityScore,ANY,missing,ANY-method}
\alias{[[,QualityScore,ANY,missing-method}
\alias{append,QualityScore,QualityScore,missing-method}
\alias{length,QualityScore-method}
\alias{width,QualityScore-method}
\alias{detail,QualityScore-method}
% NumericQuality
\alias{width,NumericQuality-method}
\alias{show,NumericQuality-method}
% MatrixQuality
\alias{[,MatrixQuality,ANY,missing-method}
\alias{[,MatrixQuality,ANY,missing,ANY-method}
\alias{[[,MatrixQuality,ANY,missing-method}
\alias{dim,MatrixQuality-method}
\alias{length,MatrixQuality-method}
\alias{narrow,MatrixQuality-method}
\alias{width,MatrixQuality-method}
\alias{append,MatrixQuality,MatrixQuality,missing-method}
% FastqQuality
\alias{width,FastqQuality-method}
\alias{show,FastqQuality-method}
\alias{alphabet,FastqQuality-method}
\alias{alphabetFrequency,FastqQuality-method}
\alias{alphabetByCycle,FastqQuality-method}
\alias{alphabetScore,FastqQuality-method}
\alias{coerce,FastqQuality,numeric-method}
\alias{coerce,FastqQuality,matrix-method}
\alias{coerce,FastqQuality,PhredQuality-method}
\alias{narrow,FastqQuality-method}
\alias{compact,FastqQuality-method}
\alias{srsort,FastqQuality-method}
\alias{srorder,FastqQuality-method}
\alias{srrank,FastqQuality-method}
\alias{srduplicated,FastqQuality-method}
% SFastqQuality
\alias{alphabetScore,SFastqQuality-method}
\alias{coerce,SFastqQuality,matrix-method}
\alias{coerce,SFastqQuality,SolexaQuality-method}

\title{Quality scores for short reads and their alignments}
\description{

  This class hierarchy represents quality scores for short
  reads. \code{QualityScore} is a virtual base class, with derived
  classes offering different ways of representing qualities. Methods
  defined on \code{QualityScore} are implemented in all derived
  classes.

}

\section{Objects from the Class}{
  Objects from the class are created using constructors (e.g.,
  \code{\link{NumericQuality}}) named after the class name.
}
\section{Extends}{
  Class \code{"\linkS4class{.ShortReadBase}"}, directly.
}
\section{Methods}{

  The following methods are defined on all \code{QualityScore} and
  derived classes:
  \describe{

    \item{[}{\code{signature(x = "QualityScore", i = "ANY", j = "missing")}}
    \item{[}{\code{signature(x = "MatrixQuality", i = "ANY", j = "missing")}:

      Subset the object, with index \code{i} indicating the reads for
      which quality scores are to be extracted. The class of the result
      is the same as the class of \code{x}. It is an error to provide
      any argument other than \code{i}.}

    \item{[[}{\code{signature(x = "QualityScore", i = "ANY", j = "ANY")}: 

      Subset the object, returning the quality score (e.g., numeric
      value) of the \code{i}th read. }

    \item{[[}{\code{signature(x = "MatrixQuality", i = "ANY", j = "ANY")}:

      Returns the vector of quality scores associated with the
      \code{i}th read.}

    \item{dim}{\code{signature(x = "MatrixQuality")}:

      The integer(2) dimension (e.g., number of reads, read width)
      represented by the quality score.}

    \item{length}{\code{signature(x = "QualityScore")}:}
    \item{length}{\code{signature(x = "MatrixQuality")}:

      The integer(1) length (e.g., number of reads) represented by the
      quality score. Note that \code{length} of \code{MatrixQuailty} is
      the number of rows of the corresponding matrix, and not the length
      of the corresponding numeric vector.}

    \item{append}{\code{signature(x = "QualityScore", values = "QualityScore", length = "missing")}:
      append \code{values} after \code{x}.}

    \item{width}{\code{signature(x = "QualityScore")}:}
    \item{width}{\code{signature(x = "NumericQuality")}:}
    \item{width}{\code{signature(x = "MatrixQuality")}:}
    \item{width}{\code{signature(x = "FastqQuality")}:

      A numeric vector with length equal to the number of quality
      scores, and value equal to the number of quality scores for each
      read. For instance, a \code{\link{FastqQuality}} will have widths
      equal to the number of nucleotides in the underlying short read.

    }

    \item{show}{\code{signature(object = "QualityScore")}:}
    \item{show}{\code{signature(object = "NumericQuality")}:}
    \item{show}{\code{signature(object = "FastqQuality")}:

      provide a brief summary of the object content.

    }

    \item{detail}{\code{signature(object = "QualityScore")}:

      provide a more detailed view of object content.

    }
  }

  The following methods are defined on specific classes:

  \describe{
    \item{alphabet}{\code{signature(x = "FastqQuality", ...)}:
      Return a character vector of valid quality characters.
    }

    \item{alphabetFrequency}{\code{signature(stringSet = "FastqQuality")}:

      Apply \code{\link[Biostrings:letterFrequency]{alphabetFrequency}}
      to quality scores, returning a matrix as described in
      \code{\link[Biostrings:letterFrequency]{alphabetFrequency}}.}

    \item{alphabetByCycle}{\code{signature(stringSet = "FastqQuality")}:

      Apply \code{\link{alphabetByCycle}} to quality scores, returning a
      matrix as described in \code{\link{alphabetByCycle}}.}

    \item{alphabetScore}{\code{signature(object = "FastqQuality")}:}
    \item{alphabetScore}{\code{signature(object = "SFastqQuality")}:

      Apply \code{\link{alphabetScore}} (i.e., summed base quality, per
      read) to \code{object}.}

    \item{coerce}{\code{signature(from = "FastqQuality", to = "numeric")}:}
    \item{coerce}{\code{signature(from = "FastqQuality", to = "matrix")}:}
    \item{coerce}{\code{signature(from = "FastqQuality", to = "PhredQuality")}:}
    \item{coerce}{\code{signature(from = "SFastqQuality", to = "matrix")}:}

    \item{coerce}{\code{signature(from = "SFastqQuality", to = "SolexaQuality")}:

      Use like \code{as(from, "matrix")}) to coerce objects of class
      \code{from} to class \code{to}, using the quality encoding implied
      by the class. When \code{to} is \dQuote{matrix}, all quality
      scores must be of the same width, and the result is a matrix of
      type \code{integer}. The result always represents the
      integer encoding of the corresponding quality string.}

    \item{narrow}{\code{signature(x = "FastqQuality", start = NA, end = NA, width = NA, use.names = TRUE)}:
      \sQuote{narrow} \code{quality} so that scores are between
      \code{start} and \code{end} bases, according to
      \code{\link[IRanges:Ranges-utils]{narrow}} in the \code{IRanges}
      package.}

    \item{compact}{\code{signature(x = "FastqQuality", ...)}:
	  reduce the space (memory) occupied by \code{x}, if
	  possible.}

    \item{srorder}{\code{signature(x = "FastqQuality")}:}
    \item{srrank}{\code{signature(x = "FastqQuality")}:}
    \item{srduplicated}{\code{signature(x = "FastqQuality")}:
    
      Apply \code{\link{srsort}}, \code{srorder}, \code{srrank}, and
      \code{srduplicated} to quality scores, returning objects as
      described on the appropriate help page.}

  }

  Integer representations of \code{SFastqQuality} and
  \code{FastqQuality} can be obtained with \code{as(x, "matrix")}.

}
\seealso{
  \code{\link{NumericQuality}} and other constructors.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\examples{
names(slot(getClass("QualityScore"), "subclasses"))
}
\keyword{classes}
