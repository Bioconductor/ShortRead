\name{readMaqMapview}
\alias{readMaqMapview}

\title{
  Read MAQ mapview text files into AlignedRead R objects
}
\description{

  The result of the MAQ command \code{maq mapview} is a text file
  summarizing alignments. This function reads one or more of these
  alignment files into a \code{\linkS4class{AlignedRead}} object.

}
\usage{
readMaqMapview(dirPath, pattern=character(0),
             sep = "\t", header = FALSE)
}

\arguments{

  \item{dirPath}{A character vector giving the directory path (relative
    or absolute) of files to be read.}

  \item{pattern}{The (\code{\link{grep}}-style) pattern describing file
    names to be read. The default (\code{character(0)}) matches all
    files in the directory.}

  \item{sep}{A length 1 character vector describing the column
    separator.}

  \item{header}{A length 1 logical vector indicating whether files
    include a header line identifying columns. If present, the header of
    the first file is used to name the returned values.}

}

\details{

  Interpretation of columns is based on the description in the MAQ
  manual, specifically

  \preformatted{
    ...each line consists of read name, chromosome, position, strand,
    insert size from the outer coorniates of a pair, paired flag,
    mapping quality, single-end mapping quality, alternative mapping
    quality, number of mismatches of the best hit, sum of qualities of
    mismatched bases of the best hit, number of 0-mismatch hits of the
    first 24bp, number of 1-mismatch hits of the first 24bp on the
    reference, length of the read, read sequence and its quality.
  }

  The read name, read sequence, and quality are read as
  \code{XStringSet} objects. Chromosome and strand are read as
  \code{factor}s.  Position is \code{numeric}, while mapping quality is
  \code{numeric}. These fields are mapped to their corresponding
  representation in \code{AlignedRead} objects.

  Number of mismatches of the best hit, sum of qualities of mismatched
  bases of the best hit, number of 0-mismatch hits of the first 24bp,
  number of 1-mismatch hits of the first 24bp are represented in the
  \code{AlignedRead} object as components of \code{alignData}.

  Remaining fields are currently ignored.

}

\value{

  A \code{\linkS4class{AlignedRead}} object.

}

\references{

  The MAQ reference manual,
  \url{http://maq.sourceforge.net/maq-manpage.shtml#5}, 3 May, 2008

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
dirPath <- system.file('extdata', 'maq', package='ShortRead')
list.files(dirPath)
countLines(dirPath)
## two files collapse into one
readMaqMapview(dirPath)
}

\keyword{IO}
