\name{ShortRead-deprecated}
\alias{uniqueFilter}

\title{Deprecated functions from the ShortRead package}

\description{

  These functions are deprecated, and will become defunct.

}

\usage{

uniqueFilter(withSread=TRUE, .name="UniqueFilter")

}

\arguments{

  \item{withSread}{A \code{logical(1)} indicating whether uniqueness
    includes the read sequence (\code{withSread=TRUE}) or is based only
    on chromosome, position, and strand (\code{withSread=FALSE})}

  \item{.name}{An optional \code{character(1)} object used to over-ride
    the name applied to default filters.}

}

\details{

  See \code{\link{srFilter}} for details of ShortRead filters.

  \code{uniqueFilter} selects elements satisfying
  \code{!srduplicated(x)} when \code{withSread=TRUE}, and
  \code{!(duplicated(chromosome(x)) & duplicated(position(x)) &
	duplicated(strand(x)))} when \code{withSread=FALSE}.

  The behavior when \code{withSread=TRUE} can be obtained with
  \code{occurrenceFilter(withSread=TRUE)}. The behavior when
  \code{withSread=FALSE} can be obtained using a custom filter

}