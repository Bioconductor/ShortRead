\name{spViewPerGene}
\alias{spViewPerGene}

\title{
Tools to visualize genomic data
}
\description{
Use \code{Snapshot}-class to visulatize a specific region of genomic data
}
\usage{
spViewPerGene(GRL, name, files, ignore.strand=FALSE, multi.lines = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GRL}{Object \code{GRangeList} containing annotation of
    genomic data. It can be generated from applying  \code{exonsBy()}
    \code{transcriptsBy()} to a \code{TrascriptDb} instance.}
  
  \item{name}{Character(1) indicating the name of element subject to be plot.}

  \item{files}{Charactor() or \code{BamFileList} specifying the file(s)
    to be visuliazed.}
  
  \item{ignore.strand}{Logical(1) indicating whether to ignore the strand of the    genomic data.}

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
txdbFile <- system.file("extdata", "sacCer2_sgdGene.sqlite",
                         package="yeastNagalakshmi")
txdb <- loadFeatures(txdbFile)
grl <- exonsBy(txdb, by="gene")
file <- system.file("extdata", "SRR002051.chrI-V.bam",
                         package="yeastNagalakshmi")
s <- spViewPerFeature(GRL=grl, name="YAL001C", files=file)

## multi-files
bfiles <- BamFileList(c(file, file))
values(bfiles) <- DataFrame(sampleGroup=factor(c("normal", "tumor")))
values(bfiles)
s <- spViewPerFeature(GRL=grl, name="YAL001C",
                      files=bfiles, multi.levels=TRUE, fac="sampleGroup")

}

