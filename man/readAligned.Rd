\name{readAligned}

\alias{readAligned}
\alias{readAligned,character-method}

\title{Read aligned reads and their quality scores into R representations}

\description{

  \code{readAligned} reads all aligned read files in a directory
  \code{dirPath} whose file name matches \code{pattern},
  returning a compact internal representation of the alignments,
  sequences, and quality scores in the files. Methods read all files into a
  single R object; a typical use is to restrict input to a single
  aligned read file.

}
\usage{

readAligned(dirPath, pattern=character(0), ...)

}

\arguments{

  \item{dirPath}{A character vector (or other object; see methods
    defined on this generic) giving the directory path (relative or
    absolute) of aligned read files to be input.}

  \item{pattern}{The (\code{\link{grep}}-style) pattern describing file
    names to be read. The default (\code{character(0)}) results in
    (attempted) input of all files in the directory.}

  \item{...}{Additional arguments, perhaps used by methods.}

}

\details{

  There is no standard aligned read file format; methods parse
  particular file types.

  The \code{readAligned,character-method} interprets file types based
  on an additional \code{type} argument. Supported types are:
  \itemize{
    \item{\code{type="SolexaExport"}}{

      This type parses \code{s_L_export.txt} files following the
      documentation in the Solexa Genome Aligment software manual,
      version 0.3.0. These files consist of the following columns;
      consult Solexa documentation for precise descriptions. If parsed,
      values can be retrieved from \code{\linkS4class{AlignedRead}} as
      follows:

      \describe{
	\item{Machine}{Ignored}
	\item{Run number}{stored in \code{alignData}}
	\item{Lane}{stored in \code{alignData}}
	\item{Tile}{stored in \code{alignData}}
	\item{X}{stored in \code{alignData}}
	\item{Y}{stored in \code{alignData}}
	\item{Index string}{Ignored}
	\item{Read number}{Ignored}
	\item{Read}{\code{sread}}
	\item{Quality}{\code{quality}}
	\item{Match chromosome}{\code{chromosome}}
	\item{Match contig}{Ignored}
	\item{Match position}{\code{position}}
	\item{Match strand}{\code{strand}}
	\item{Match description}{Ignored}
	\item{Single-read alignment score}{\code{alignQuality}}
	\item{Paired-read alignment score}{Ignored}
	\item{Partner chromosome}{Ignored}
	\item{Partner contig}{Ignored}
	\item{Partner offset}{Ignored}
	\item{Partner strand}{Ignored}
	\item{Filtering}{\code{alignData}}
      }

      Paired read columns are not interpretted.  The resulting
      \code{\linkS4class{AlignedRead}} object does \emph{not} contain a
      meaningful \code{id}; instead, use information from
      \code{alignData} to identify reads.

      Different interfaces to reading alignment files are described in
      \code{\linkS4class{SolexaPath}} and \code{\linkS4class{SolexaSet}}.

    }

    \item{\code{type="MAQMapview"}}{

      This type parses alignment files created by MAQ.  Interpretation
      of columns is based on the description in the MAQ manual,
      specifically

      \preformatted{

        ...each line consists of read name, chromosome, position,
        strand, insert size from the outer coorniates of a pair,
        paired flag, mapping quality, single-end mapping quality,
        alternative mapping quality, number of mismatches of the
        best hit, sum of qualities of mismatched bases of the best
        hit, number of 0-mismatch hits of the first 24bp, number
        of 1-mismatch hits of the first 24bp on the reference,
        length of the read, read sequence and its quality.

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
  }
}

\value{

  A single R object (e.g., \code{\linkS4class{AlignedRead}}) containing
  alignments, sequences and qualities of all files in \code{dirPath}
  matching \code{pattern}. There is no guarantee of order in which files
  are read.

}

\seealso{

  A \code{\linkS4class{AlignedRead}} object.

  The MAQ reference manual,
  \url{http://maq.sourceforge.net/maq-manpage.shtml#5}, 3 May, 2008

}

\author{Martin Morgan}

\examples{
dirPath <- system.file('extdata', 'maq', package='ShortRead')
list.files(dirPath)
## First line
readLines(list.files(dirPath, full.names=TRUE)[[1]], 1)
countLines(dirPath)
## two files collapse into one
readAligned(dirPath, type="MAQMapview")
}
\keyword{manip}
