\name{countLines}

\alias{countLines}

\title{Count lines in all (text) files in a directory whose file name
  matches a pattern}

\description{

  \code{countLines} visits all files in a directory path \code{dirPath}
  whose base (i.e., file) name matches \code{pattern}. Lines in the file
  are counted as the number of new line characters.

}
\usage{

countLines(dirPath, pattern=character(0), ..., useFullName=FALSE)

}

\arguments{

  \item{dirPath}{A character vector (or other object; see methods
    defined on this generic) giving the directory path (relative or
    absolute) of files whose lines are to be counted.}

  \item{pattern}{The (\code{\link{grep}}-style) pattern describing files
    whose lines are to be counted. The default (\code{character(0)})
    results in line counts for all files in the directory.}

  \item{...}{Additional arguments, passed internally to list.files. See
    \code{\link{list.files}}.}

  \item{useFullName}{A \code{logical(1)} indicating whether elements of
    the returned vector should be named with the base (file) name
    (default; \code{useFullName=FALSE}) or the full path name
    (\code{useFullName=TRUE}).}

}
\value{

  A named integer vector of line counts. Names are paths to the files
  whose lines have been counted, excluding \code{dirPath}.

}
\author{Martin Morgan}
\examples{
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
countLines(analysisPath(sp))
countLines(experimentPath(sp), recursive=TRUE)
countLines(experimentPath(sp), recursive=TRUE, useFullName=TRUE)
}
\keyword{manip}
