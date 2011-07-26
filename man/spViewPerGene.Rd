\name{spViewPerGene}
\alias{spViewPerGene}

\title{
Tools to visualize genomic data
}
\description{
Use \code{Snapshot}-class to visulatize a specific region of genomic data
}
\usage{
spViewPerGene(GRL, name, files, multi.lines = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GRL}{Object \code{GRangeList} containing annotation of
    genomic data. It can be generated from applying  \code{exonsBy()}
    \code{transcriptsBy()} to a \code{TrascriptDb} instance.}
  
  \item{name}{Character(1) indicating the name of element subject to be plot.}

  \item{files}{Charactor() or \code{BamFileList} specifying the file(s)
    to be visuliazed.}

  \item{multi.lines}{Logical(1) inidicating whether to see multiple
    files on different panels. If \code{FALSE}, }
  \item{\dots}{Argument used for createing a \code{Snapshot} object.}
}

\value{A \code{Snapshot} instance}

\author{Chao-Jen Wong \email{cwon2@fhcrc.org}}

\seealso{
\code{\link{Snapshot}}
}
\examples{
library(GenomicFeatures)

}

