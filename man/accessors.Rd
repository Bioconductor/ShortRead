\name{accessors}

\alias{accessors}
% SRVector
\alias{vclass}
% ShortRead / ShortReadQ
\alias{sread}
\alias{id}
\alias{quality}
% AlignedRead
\alias{chromosome}
\alias{position}
\alias{strand}
\alias{alignQuality}
\alias{alignData}
% Solexa
\alias{experimentPath}
\alias{dataPath}
\alias{scanPath}
\alias{imageAnalysisPath}
\alias{baseCallPath}
\alias{analysisPath}
% SolexaSet
\alias{solexaPath}
\alias{laneDescription}
\alias{laneNames}

\title{Accessors for ShortRead classes}

\description{

  These functions are `accessors' (to get and set values) for objects in
  the \pkg{ShortRead} package.

}
\usage{

## SRVector
vclass(object, ...)
## ShortRead / ShortReadQ
sread(object, ...)
id(object, ...)
quality(object, ...)
## AlignedRead
chromosome(object, ...)
position(object, ...)
strand(object, ...)
alignQuality(object, ...)
alignData(object, ...)
## Solexa
experimentPath(object, ...)
dataPath(object, ...)
scanPath(object, ...)
imageAnalysisPath(object, ...)
baseCallPath(object, ...)
analysisPath(object, ...)
## SolexaSet
solexaPath(object, ...)
laneDescription(object, ...)
laneNames(object, ...)

}

\arguments{

  \item{object}{An object derived from class \code{ShortRead}. See
    help pages for individual objects, e.g.,
    \code{\linkS4class{ShortReadQ}}. The default is to extract the
    contents of a slot of the corresponding name (e.g., slot
    \code{sread}) from \code{object}.}

  \item{...}{Additional arguments passed to the accessor. The default
    definitions do not make use of additional arguments.}

}
\value{

  Usually, the value of the corresponding slot, or other simple content
  described on the help page of \code{object}.

}
\author{Martin Morgan}
\examples{
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
experimentPath(sp)
basename(analysisPath(sp))
}
\keyword{manip}
