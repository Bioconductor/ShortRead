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
\alias{append,QualityScore,QualityScore-method}
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
\alias{append,MatrixQuality,MatrixQuality-method}
% FastqQuality
\alias{width,FastqQuality-method}
\alias{show,FastqQuality-method}
\alias{alphabet,FastqQuality-method}
\alias{encoding,FastqQuality-method}
\alias{alphabetFrequency,FastqQuality-method}
\alias{alphabetByCycle,FastqQuality-method}
\alias{alphabetScore,FastqQuality-method}
\alias{coerce,FastqQuality,numeric-method}
\alias{coerce,FastqQuality,matrix-method}
\alias{coerce,FastqQuality,PhredQuality-method}
\alias{narrow,FastqQuality-method}
\alias{trimTailw,FastqQuality-method}
\alias{trimTails,FastqQuality-method}
\alias{srsort,FastqQuality-method}
\alias{srorder,FastqQuality-method}
\alias{srrank,FastqQuality-method}
\alias{srduplicated,FastqQuality-method}
% SFastqQuality
\alias{encoding,SFastqQuality-method}
\alias{alphabetScore,SFastqQuality-method}
\alias{coerce,SFastqQuality,matrix-method}
\alias{coerce,SFastqQuality,SolexaQuality-method}

%% Biostrings::PhredQuality
\alias{alphabetScore,PhredQuality-method}


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

  Defined classes are as follows:
  \describe{

    \item{QualityScore}{Virtual base class; instances cannot be
      instantiated.}

    \item{NumericQuality}{A single numeric vector, where values
      represent quality scores on an arbitrary scale.}

    \item{IntegerQuality}{A integer numeric vector, where values
      represent quality scores on an arbitrary scale.}

    \item{MatrixQuality}{A rectangular matrix of quality scores, with
      rows representing reads and columns cycles. The content and
      interpretation of row and column entries is arbitrary; the
      rectangular nature implies quality scores from equal-length
      reads.}

    \item{FastqQuality}{\sQuote{fastq} encoded quality scores stored in
      a \code{BStringSet} instance. Base qualities of a single read are
      represented as an ASCII character string. The integer-valued
      quality score of a single base is encoded as its ASCII equivalent
      plus 33. The precise definition of the integer-valued quality
      score is unspecified, but is usually a Phred score; the meaning
      can be determined from the source of the quality scores. Multiple
      reads are stored as a \code{BStringSet}, and so can be of varying
      lengths.}

    \item{SolexaQuality}{As with \code{FastqQuality}, but with integer
      qualities encoded as ASCII equivalent plus 64.}
    
  }

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

    \item{append}{\code{signature(x = "QualityScore", values = "QualityScore")}:
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

    \item{detail}{\code{signature(x = "QualityScore")}:

      provide a more detailed view of object content.

    }
  }

  The following methods are defined on specific classes:

  \describe{
    \item{alphabet}{\code{signature(x = "FastqQuality", ...)}:
      Return a character vector of valid quality characters.
    }

    \item{encoding}{\code{signature(x = "FastqQuality", ...)},
      \code{signature(x = "SFastqQuality", ...)}:
      Returns a named character vector of integer encodings.
    }

    \item{alphabetFrequency}{\code{signature(stringSet = "FastqQuality")}:

      Apply \code{\link[Biostrings:letterFrequency]{alphabetFrequency}}
      to quality scores, returning a matrix as described in
      \code{\link[Biostrings:letterFrequency]{alphabetFrequency}}.}

    \item{alphabetByCycle}{\code{signature(stringSet = "FastqQuality")}:

      Apply \code{\link{alphabetByCycle}} to quality scores, returning a
      matrix as described in \code{\link{alphabetByCycle}}.}

    \item{alphabetScore}{\code{signature(object = "FastqQuality")}:}
    \item{alphabetScore}{\code{signature(object = "SFastqQuality")}:}
    \item{alphabetScore}{\code{signature(object = "PhredQuality")}:

      Apply \code{\link{alphabetScore}} (i.e., summed base quality, per
      read) to \code{object}.}

    \item{coerce}{\code{signature(from = "FastqQuality", to = "numeric")}:}
    \item{coerce}{\code{signature(from = "FastqQuality", to = "matrix")}:}
    \item{coerce}{\code{signature(from = "FastqQuality", to = "PhredQuality")}:}
    \item{coerce}{\code{signature(from = "SFastqQuality", to = "matrix")}:}

    \item{coerce}{\code{signature(from = "SFastqQuality", to = "SolexaQuality")}:

      Use \code{as(from, "matrix")}) and similar to coerce objects of
      class \code{from} to class \code{to}, using the quality encoding
      implied by the class. When \code{to} is \dQuote{matrix}, the
      result is a matrix of type \code{integer} with number of columns
      equal to the maximum width of \code{from}; elements \code{i, j}
      with \code{j > width(from)[i]} have value \code{NA_integer_}. The
      result always represents the integer encoding of the corresponding
      quality string.}

    \item{narrow}{\code{signature(x = "FastqQuality", start = NA, end = NA, width = NA, use.names = TRUE)}:
      \sQuote{narrow} \code{quality} so that scores are between
      \code{start} and \code{end} bases, according to
      \code{\link[IRanges:intra-range-methods]{narrow}} in the \code{IRanges}
      package.}

    \item{trimTailw}{\code{signature(object="FastqQuality", k="integer",
	a="character", halfwidth="integer", ..., ranges=FALSE)}:
      trim trailing nucleotides when a window of width 2 * halfwidth + 1
      contains \code{k} or more quality scores falling at or below
      \code{a}.}

    \item{trimTails}{\code{signature(object="FastqQuality", k="integer",
	a="character", successive=FALSE, ..., ranges=FALSE)}: trim
      trailing scores if \code{k} scores fall below the quality encoded
      by \code{a}. If \code{successive=FALSE}, the k'th failing score
      and all subseqent scores are trimmed. If \code{successive=TRUE},
      failing scores must occur successively; the sequence is trimmed
      from the first of the successive failing score.}

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
encoding(FastqQuality())
encoding(SFastqQuality())
}
\keyword{classes}
