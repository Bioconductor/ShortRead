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

\title{Class "SnapshotFunction"}

\description{
A class to store custom reader and viewer functions for the
\code{\link{Snapshot}} class.
}

\usage{
SnapshotFunction(reader, viewer, limits, ...)
}

\arguments{
  \item{reader}{A \code{\link{function}} for reading data from BAM files and
      returning a \code{\link{data.frame}.}}
  \item{viewer}{A \code{\link{function}} for visulazation that returns an
      \code{\link{SpTrellis}} object.}
  \item{limits}{Limits of the ranges to be visualized.}
  \item{...}{Additional arguments to be passed to the SnapshotFunction function.}
}

\section{Fields}{
  \describe{
    \item{\code{reader}:}{Object of class \code{\link{function}} for reading data from BAM files and
      returning a \code{\link{data.frame}}.}
    \item{\code{viewer}:}{Object of class \code{\link{function}} for visulazation that returns an
      \code{\link{SpTrellis}} object.}
    \item{\code{...}:}{Object of class \code{GRanges}.}
    \item{\code{limits}:}{Object of class \code{integer} for the limits of
        ranges to be visuzlized.}
  }
}

\author{Chao-Jen \url{cwon2@fhcrc.org}}

\seealso{\code{\link{Snapshot}}}

\examples{showClass("SnapshotFunction")}

\keyword{classes}

