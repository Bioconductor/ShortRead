\name{srdistance}
\alias{srdistance}
%
\alias{srdistance,DNAStringSet,character-method}
\alias{srdistance,DNAStringSet,DNAString-method}
\alias{srdistance,DNAStringSet,DNAStringSet-method}

\title{Edit distances between reads and a small number of short references}

\description{

  \code{srdistance} calculates the edit distance from each read in
  \code{pattern} to each read in \code{subject}. The underlying
  algorithm \code{\link[Biobase]{pairwiseAlignment}} is only efficient when both
  reads are short, and when the number of \code{subject} reads is small.

}

\usage{
srdistance(pattern, subject, ...)
}

\arguments{
  \item{pattern}{An object of class \code{DNAStringSet} containing reads
    whose edit distance is desired.}

  \item{subject}{A short \code{character} vector, \code{DNAString} or
    (small) \code{DNAStringSet} to serve as reference.}

  \item{\dots}{additional arguments, forward to \code{srapply}.}
}

\details{

  The underlying algorithm performs pairwise alignment from each read in
  \code{pattern} to each sequence in \code{subject}. The return value is
  a list of numeric vectors of distances, one list element for each
  sequence in \code{subject}. The vector in each list element contains
  for each read in \code{pattern} the edit distance from the read to the
  corresponding subject. The weight matrix and gap penalties used to
  calculate the distance are structured to weight base substitutions and
  single base insert/deletions equally. Edit distance between known and
  ambiguous (e.g., N) nucleotides, or between ambiguous nucleotides, are
  weighted as though each possible nucleotide in the ambiguity were
  equally likely.

}

\value{

  A list of length equal to that of \code{subject}. Each element is a
  numeric vector equal to the length of \code{pattern}, with values
  corresponding to the minimum distance between between the
  corresponding pattern and subject sequences.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{\code{\link[Biobase]{pairwiseAlignment}}}

\examples{
sp <- SolexaPath(system.file("extdata", package="ShortRead"))
aln <- readAligned(sp, "s_2_export.txt")
polyA <- polyn("A", 35)
polyT <- polyn("T", 35)

d1 <- srdistance(clean(sread(aln)), polyA)
d2 <- srdistance(sread(aln), polyA)
d3 <- srdistance(sread(aln), c(polyA, polyT))
}

\keyword{manip}
