\name{AlignedRead}
\alias{AlignedRead}

\title{Construct objects of class "AlignedRead"}
\description{

  This function constructs objects of
  \code{\linkS4class{AlignedRead}}. It will often be more convenient to
  create \code{AlignedRead} objects using parsers such as
  \code{\link{readMaqMapview}}.

}
\usage{
AlignedRead(sread, id, quality, chromosome, position, strand,
            alignQuality,
            alignData = AlignedDataFrame(nrow = length(sread)))
}

\arguments{
  \item{sread}{An object of class \code{DNAStringSet}, containing the
    DNA sequences of the short reads.}

  \item{id}{An object of class \code{BStringSet}, containing the
    identifiers of the short reads. This object is the same length as
    \code{sread}.}

  \item{quality}{An object of class \code{BStringSet}, containing the
    ASCII-encoded quality scores of the short reads. This object is the
    same length as \code{sread}.}

  \item{chromosome}{A \code{factor} describing the chromosome to which
    each short read aligns.}

  \item{position}{A \code{integer} vector describing the (base pair)
    position at which each short read begins its alignment.}

  \item{strand}{A \code{factor} describing the strand to which the short
    read aligns.}

  \item{alignQuality}{A \code{numeric} vector describing the alignment
    quality.}

  \item{alignData}{An \code{AlignedDataFrame} with number of rows equal
    to the length of \code{sread}, containing additional information
    about alignments.}
}
\value{
  An object of class \code{\linkS4class{AlignedRead}}.
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{\code{\linkS4class{AlignedRead}}.}
\keyword{manip}
