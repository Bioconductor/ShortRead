\name{pileup}
\alias{pileup}

\title{Calculate a pile-up representation of short-read mappings}

\description{
Given short read mappings or similar data, this function calculates a pile-up, i.e. 
representing the reference sequence (that is, typically, one of the chromosome), such
that its length is the number of base pairs of the reference sequence, and each integer 
is the number of reads (or fragments, see below) mapped to the corresponding basepair. 
}

\usage{
pileup( start, fraglength, chrlength, 
   dir = factor( "+", levels=c("-","+") ), 
   readlength = fraglength,
   offset = 1 )
}

\arguments{

\item{start}{A vector with the start positions of each read on the reference sequence. All
reads must correspond to the same reference sequence.}

\item{fraglength}{A vector of the same length as 'start' with the lengths of all the fragments.
Alternatively, a single integer, specifying one constant length to assume for all tags.}

\item{chrlength}{The length of the reference sequence. You may use the function
\code{\link{readBfaToc}} to extract this information from the .bfa file.}

\item{dir}{A factor with level "-" and "+" of the same length as 'start', specifying whether the
fragment extends to the right (towards higher index values, '+') or to the left (towards lower
index values, '-') beyond the read. See below for more explanation.}

\item{readlength}{The length of the reads, either as a vector of the same length as 'start' or as
a single number. This parameter makes sense only if 'dir' is used, too. If not specified, read
lengths and fragment lengths are taken to be the same.}

\item{offset}{The index of the first base pair in the result vector. The default is 1, i.e.
assumes that the 'start' positions are in 1-based chromosome coordinates.}

}

\value{an integer vector of length 'chrlength', each element counting how many fragments map to
this basepair.}

\note{

1. This function is not (yet) suitable for paired-end reads.


2. If the arguments 'dir' and 'readlength' are not used, the fragments are assumed to start at the
positions given in 'start' and extend to the right by the number of basepairs given in
fraglength. If 'dir' and 'readlength' are supplied then the interval starting at 'start' and
extending to the right by the number of base pairs given in 'readlength' marks the position of
the read, which is one end of the fragment. If 'dir' ist '+', it is taken as the left end and 
the fragment will be extended to the right to have the total length given by 'fraglength'. If
'dir' is '-', the end is taken as the right end and is extended to the left. Note that in the
latter case, the 'start' position does mark the border between read and rest of fragment, not an
actual 'end' of the fragment. If you are confused now, look at the examples below.

3. Sorry for the inconsequent use of 'width' and 'length' in a seemingly interchangeable fashion.
}

\examples{\dontrun{ 

Example 1: Assuming that 'lane' is a \link[=AlignedRead-class]{AlignedRead} object containing aligned reads froma
Solexa lane, you may get a pile-up representation of chromosome 13 as follows

chr13length <- 114142980   # the length of human chromosme 13
pu <- pileup( position(lane)[chromosome(lane)=="13"], width(lane), chr13length )

Example 2: Even though the width of the reads (as repored by \code{width(lane)}) is only 24,
these 24 bp are just one end of a longer fragment. Assuming that all fragments have been
sonicated to about the same length, say 150 bp, we may get a better pile-up representation by:

pu2 <- pileup( position(lane)[chromosome(lane)=="13"], 150, chr13length,
strand(lane)[chromosome(lane)=="13"], width(lane) )

}}

\author{Simon Anders, EMBL-EBI, \email{sanders@fs.tum.de}}
