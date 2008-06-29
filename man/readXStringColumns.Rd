\name{readXStringColumns}
\alias{readXStringColumns}

\title{
  Read one or more columns into XStringSet (e.g., DNAStringSet) objects
}

\description{

  This function allows short read data components such as DNA sequence,
  quality scores, and read names to be read in to \code{XStringSet}
  (e.g., \code{DNAStringSet}, \code{BStringSet}) objects. One or several
  files of identical layout can be specified.

}
\usage{
readXStringColumns(dirPath, pattern=character(0),
                   colClasses=list(NULL), sep = "\t",
                   header = FALSE, comment.char="#")
}
\arguments{

  \item{dirPath}{A character vector giving the directory path (relative
    or absolute) of files to be read.}

  \item{pattern}{The (\code{\link{grep}}-style) pattern describing file
    names to be read. The default (\code{character(0)}) reads all files
    in \code{dirPath}. All files are expected to have identical numbers
    of columns.}

  \item{colClasses}{A list of length equal to the number of columns in
    a file. Columns with corresponding \code{coClasses} equal to
    \code{NULL} are ignored. Other entries in \code{colClasses} are
    expected to be character strings describing the base class for the
    \code{XStringSet}. For instance a column of DNA sequences would be
    specified as \code{"DNAString"}. The column would be parsed into a
    \code{DNAStringSet} object.}

  \item{sep}{A length 1 character vector describing the column
    separator.}

  \item{header}{A length 1 logical vector indicating whether files
    include a header line identifying columns. If present, the header of
    the first file is used to name the returned values.}

  \item{comment.char}{A length 1 character vector, with a single
    character that, when appearing at the start of a line, indicates
    that the entire line should be ignored. Currently there is no way to
    use comment characters in other than the first position of a line.}
}

\value{

  A list, with each element containing an \code{XStringSet} object of
  the type corresponding to the non-NULL elements of \code{colClasses}. 

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
## valid character strings for colClasses
names(slot(getClass("XString"), "subclasses"))

dirPath <- system.file('extdata', 'maq', package='ShortRead')

colClasses <- rep(list(NULL), 16)
colClasses[c(1, 15, 16)] <- c("BString", "DNAString", "BString")

## read one file
readXStringColumns(dirPath, "out.aln.1.txt", colClasses=colClasses)

## read all files into a single object for each column
res <- readXStringColumns(dirPath, colClasses=colClasses)
}

\keyword{IO}
