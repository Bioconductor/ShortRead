\name{SnapshotFunction-class}

\Rdversion{1.1}
\docType{class}

\alias{SnapshotFunction-class}
\alias{SnapshotFunction}
\alias{show,SnapshotFunction-method}
\alias{SnapshotFunctionList-class}
\alias{SnapshotFunctionList}
\alias{SnapshotFunctionList,Snapshot-method}
\alias{SnapshotFunctionList,ANY-method}
\alias{reader}
\alias{viewer}
\alias{limits}

\title{Class "SnapshotFunction"}

\description{

  A class to store custom reader and viewer functions for the
  \code{\link{Snapshot}} class.

}

\usage{
SnapshotFunction(reader, viewer, limits, ...)
reader(x, ...)
viewer(x, ...)
limits(x, ...)
}

\arguments{

  \item{reader}{A function for reading data. The function must take a
    single argument (a \code{\link{Snapshot}} instance) and return a
    \code{data.frame} summarizing the file.}

  \item{viewer}{A function for visualizing the data. The function must
    accept the \code{data.frame} created by \code{reader}, and return
    an \code{\link{SpTrellis}} object representing the view.}

  \item{limits}{An integer(2) indicating the minimum and maximum number
    of nucleotides the \code{SnapshotFunction} is intended to
    visualize. For instance, a \sQuote{fine-scale} viewer displaying a
    pileup might be appropriate at between 1000 and 50000 nucleotides.}

  \item{x}{An instance of \code{SnapshotFunction}}

  \item{...}{Additional arguments, currently unused.}

}

\section{Fields}{
  \describe{

    \item{\code{reader}:}{Object of class \code{\link{function}} for
      reading data from BAM files and returning a
      \code{\link{data.frame}}.}

    \item{\code{viewer}:}{Object of class \code{\link{function}} for
      visulazation that returns an \code{\link{SpTrellis}} object.}

    \item{\code{limits}:}{Object of class \code{integer} for the limits
      of ranges to be visuzlized.}
  }
}

\author{Chao-Jen \url{cwon2@fhcrc.org}}

\seealso{\code{\link{Snapshot}}}

\examples{
## internally defined function
reader(ShortRead:::.fine_coverage)
viewer(ShortRead:::.fine_coverage)
limits(ShortRead:::.fine_coverage)
}

\keyword{classes}

