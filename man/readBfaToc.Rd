\name{readBfaToc}
\alias{readBfaToc}

\title{(Legacy) Get a list of the sequences in a Maq .bfa file}

\description{
As \code{\link{coverage}} needs to know the lengths of the reference sequences,
this function is provided which extracts this information from
a .bfa file (Maq's "binary FASTA" format).
}

\usage{

readBfaToc( bfafile )
}

\arguments{
\item{bfafile}{The file name of the .bfa file.}
}

\value{An integer vector with one element per
reference sequence found in the .bfa file, each vector element named with the
sequence name and having the sequence length as value.}

\author{Simon Anders, EMBL-EBI, \email{sanders@fs.tum.de}

(Note: The C code for this function incorporates code from Li Heng's MAQ
software, (c) Li Heng and released by him under GPL 2.}
