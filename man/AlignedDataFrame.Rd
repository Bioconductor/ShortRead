\name{AlignedDataFrame}
\alias{AlignedDataFrame}

\title{AlignedDataFrame constructor}
\description{

  Construct an \code{AlignedDataFrame} from a data frame and its
  metadata

}

\usage{
AlignedDataFrame(data, metadata, nrow = nrow(data))
}

\arguments{

  \item{data}{A data frame containing alignment information.}

  \item{metadata}{A data frame describing the columns of \code{data},
    and with number of rows of \code{metadata} corresponding to number
    of columns of \code{data}. . The data frame must contain a column
    \code{labelDescription} providing a verbose describtion of each
    column of \code{data}.}

  \item{nrow}{An optional argument, to be used when \code{data} is not
    provided, to construct an AlignedDataFrame with the specified number
    of rows.}
}
\value{

  An object of \code{\linkS4class{AlignedDataFrame}}.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\keyword{manip}
