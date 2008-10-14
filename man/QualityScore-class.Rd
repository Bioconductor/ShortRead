\name{QualityScore-class}
\docType{class}
\alias{QualityScore-class}
\alias{NumericQuality-class}
\alias{IntegerQuality-class}
\alias{MatrixQuality-class}
\alias{FastqQuality-class}
\alias{SFastqQuality-class}
% QualityScore
\alias{[,QualityScore,ANY,missing,ANY-method}
\alias{[[,QualityScore,ANY,missing-method}
\alias{length,QualityScore-method}
\alias{width,QualityScore-method}
\alias{detail,QualityScore-method}
% NumericQuality
\alias{width,NumericQuality-method}
\alias{show,NumericQuality-method}
% MatrixQuality
\alias{[,MatrixQuality,ANY,missing,ANY-method}
\alias{[[,MatrixQuality,ANY,missing-method}
\alias{length,MatrixQuality-method}
\alias{width,MatrixQuality-method}
% FastqQuality
\alias{width,FastqQuality-method}
\alias{show,FastqQuality-method}
\alias{alphabet,FastqQuality-method}
\alias{alphabetFrequency,FastqQuality-method}
\alias{alphabetByCycle,FastqQuality-method}
\alias{alphabetScore,FastqQuality-method}
\alias{srsort,FastqQuality-method}
\alias{srorder,FastqQuality-method}
\alias{srrank,FastqQuality-method}
\alias{srduplicated,FastqQuality-method}
% SFastqQuality
\alias{alphabetScore,SFastqQuality-method}

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

    \item{length}{\code{signature(x = "QualityScore")}:}
    \item{length}{\code{signature(x = "MatrixQuality")}:

      The numeric(1) length (e.g., number of reads) represented by the
      quality score.

    }

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
    \item{alphabet}{\code{signature(x = "FastqQuality")}:
      Return a character vector of valid quality characters.
    }

    \item{alphabetFrequency}{\code{signature(stringSet = "FastqQuality")}:

      Apply \code{\link[Biostrings]{alphabetFrequency}} to quality
      scores, returning a matrix as described in
      \code{\link[Biostrings]{alphabetFrequency}}.}

    \item{alphabetByCycle}{\code{signature(stringSet = "FastqQuality")}:

      Apply \code{\link{alphabetByCycle}} to quality scores, returning a
      matrix as described in \code{\link{alphabetByCycle}}.}

    \item{alphabetScore}{\code{signature(object = "FastqQuality")}:}
    \item{alphabetScore}{\code{signature(object = "SFastqQuality")}:

      Apply \code{\link{alphabetScore}} (i.e., summed base quality, per
      read) to \code{object}.}

    \item{srsort}{\code{signature(x = "FastqQuality")}:}
    \item{srorder}{\code{signature(x = "FastqQuality")}:}
    \item{srrank}{\code{signature(x = "FastqQuality")}:}
    \item{srduplicated}{\code{signature(x = "FastqQuality")}:
    
      Apply \code{\link{srsort}}, \code{srorder}, \code{srrank}, and
      \code{srduplicated} to quality scores, returning objects as
      described on the appropriate help page.}

  }
}
\seealso{
  \code{\link{NumericQuality}} and other constructors.
}
\author{Martin Morgan <mtmmorgan@fhcrc.org>}
\examples{
names(slot(getClass("QualityScore"), "subclasses"))
}
\keyword{classes}
