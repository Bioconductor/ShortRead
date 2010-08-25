\name{readFasta}

\alias{readFasta}
\alias{readFasta,character-method}
\alias{writeFasta}

\title{Read and write FASTA files to or from ShortRead objects}

\description{

  \code{readFasta} reads all FASTa-formated files in a directory
  \code{dirPath} whose file name matches pattern \code{pattern},
  returning a compact internal representation of the sequences and
  quality scores in the files. Methods read all files into a single R
  object; a typical use is to restrict input to a single FASTQ file.

  \code{writeFasta} writes an object to a single \code{file}, using
  \code{mode="w"} (the default) to create a new file or \code{mode="a"}
  append to an existing file. Attempting to write to an existing file
  with \code{mode="w"} results in an error.

}
\usage{
readFasta(dirPath, pattern = character(0), ...,
    nrec=-1L, skip=0L)
\S4method{readFasta}{character}(dirPath, pattern = character(0), ...,
    nrec=-1L, skip=0L)
writeFasta(object, file, ...)
}

\arguments{

  \item{dirPath}{A character vector giving the directory path (relative
    or absolute) or single file name of FASTA files to be read.}

  \item{pattern}{The (\code{\link{grep}}-style) pattern describing file
    names to be read. The default (\code{character(0)}) results in
    (attempted) input of all files in the directory.}

  \item{object}{An object to be output in \code{fasta} format.}

  \item{file}{A length 1 character vector providing a path to a file to
    the object is to be written to.}

  \item{...}{Additional arguments, used by \code{\link{writeFASTA}} or
    methods.}
  
  \item{nrec}{See \code{?read.DNAStringSet}.}

  \item{skip}{See \code{?read.DNAStringSet}.}
}

\value{

  \code{readFasta} returns a \code{\linkS4class{DNAStringSet}}.
  containing sequences and qualities contained in all files in
  \code{dirPath} matching \code{pattern}. There is no guarantee of order
  in which files are read.

  \code{writeFasta} is invoked primarily for its side effect, creating
  or appending to file \code{file}. The function returns, invisibly, the
  length of \code{object}, and hence the number of records written.

}

\author{Martin Morgan}

\examples{
showMethods("readFasta")

f1 <- system.file("extdata", "someORF.fa", package="Biostrings")

rfa <- readFasta(f1)
sread(rfa)
id(rfa)

file <- tempfile()
writeFasta(rfa, file)
readLines(file, 8)
}
\keyword{manip}
