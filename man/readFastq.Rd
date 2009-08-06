\name{readFastq}

\alias{readFastq}
\alias{writeFastq}
\alias{readFastq,character-method}

\title{Read and write FASTQ-formatted files}

\description{

  \code{readFastq} reads all FASTQ-formated files in a directory
  \code{dirPath} whose file name matches pattern \code{pattern},
  returning a compact internal representation of the sequences and
  quality scores in the files. Methods read all files into a single R
  object; a typical use is to restrict input to a single FASTQ file.

  \code{writeFastq} writes an object to a single \code{file}, using
  \code{mode="w"} (the default) to create a new file or \code{mode="a"}
  append to an existing file. Attempting to write to an existing file
  with \code{mode="w"} results in an error.

}
\usage{
readFastq(dirPath, pattern=character(0), ...)
\S4method{readFastq}{character}(dirPath, pattern=character(0),
    ..., qualityType=c("FastqQuality", "SFastqQuality"),
    filter=srFilter(), withIds=TRUE)
writeFastq(object, file, mode="w", ...)
}

\arguments{

  \item{dirPath}{A character vector (or other object; see methods
    defined on this generic) giving the directory path (relative or
    absolute) of FASTQ files to be read.}

  \item{pattern}{The (\code{\link{grep}}-style) pattern describing file
    names to be read. The default (\code{character(0)}) results in line
    (attempted) input of all files in the directory.}

  \item{object}{An object to be output in \code{fastq} format. For
    methods, use \code{showMethods(object, where=getNamespace("ShortRead"))}.}

  \item{file}{A length 1 character vector providing a path to a file to
    the object is to be written to.}

  \item{mode}{A length 1 character vector equal to either \sQuote{w} or
    \sQuote{a} to write to a new file or append to an existing file,
    respectively.}

  \item{...}{Additional arguments, perhaps used by methods.}

  \item{qualityType}{Representation to be used for quality scores, must
	be one of \code{FastqQuality} (Phred-like encoding),
	\code{SFastqQuality} (Illumina encoding).}

  \item{filter}{An object of class \code{\link{srFilter}}, used to
	filter objects of class \code{\linkS4class{ShortReadQ}} at input.}

  \item{withIds}{\code{logical(1)} indicating whether identifiers should
	be read from the fastq file.}

}

\details{

  The fastq format is not quite precisely defined. The basic definition
  used here parses the following four lines as a single record:

  \preformatted{
    @HWI-EAS88_1_1_1_1001_499
    GGACTTTGTAGGATACCCTCGCTTTCCTTCTCCTGT
    +HWI-EAS88_1_1_1_1001_499
    ]]]]]]]]]]]]Y]Y]]]]]]]]]]]]VCHVMPLAS
  }

  The first and third lines are identifiers preceded by a specific
  character (the identifiers are identical, in the case of Solexa). The
  second line is an upper-case sequence of nucleotides. The parser
  recognizes IUPAC-standard alphabet (hence ambiguous nucleotides),
  coercing \code{.} to \code{-} to represent missing values. The final
  line is an ASCII-encoded representation of quality scores, with one
  ASCII character per nucleotide.

  The encoding implicit in Solexa-derived fastq files is that each
  character code corresponds to a score equal to the ASCII character
  value minus 64 (e.g., ASCII \code{@} is decimal 64, and corresponds to
  a Solexa quality score of 0). This is different from BioPerl, for
  instance, which recovers quality scores by subtracting 33 from the
  ASCII character value (so that, for instance, \code{!}, with decimal
  value 33, encodes value 0).

  The BioPerl description of fastq asserts that the first character of
  line 4 is a \code{!}, but the current parser does not support this
  convention.

  \code{writeFastq} creates files following the specification outlined
  above, using the IUPAC-standard alphabet (hence, sequences containing
  \sQuote{.} when read will be represented by \sQuote{-} when written).

}

\value{

  \code{readFastq} returns a single R object (e.g.,
  \code{\linkS4class{ShortReadQ}}) containing sequences and qualities
  contained in all files in \code{dirPath} matching
  \code{pattern}. There is no guarantee of order in which files are
  read.

  \code{writeFastq} is invoked primarily for its side effect, creating
  or appending to file \code{file}. The function returns, invisibly, the
  length of \code{object}, and hence the number of records written.

}

\seealso{

  The IUPAC alphabet in Biostrings.

  \url{http://www.bioperl.org/wiki/FASTQ_sequence_format} for the
  BioPerl definition of fastq.

  Solexa documentation `Data analysis - documentation : Pipeline output
  and visualisation'.

}

\author{Martin Morgan}

\examples{
showMethods("readFastq")

sp <- SolexaPath(system.file('extdata', package='ShortRead'))
rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
sread(rfq)
id(rfq)
quality(rfq)

## SolexaPath method 'knows' where FASTQ files are placed
rfq1 <- readFastq(sp, pattern="s_1_sequence.txt")
rfq1

file <- tempfile()
writeFastq(rfq, file)
readLines(file, 8)
}
\keyword{manip}
