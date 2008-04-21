\name{detail}

\alias{detail}

\title{Show (display) detailed object content}

\description{

  This is a variant of \code{\link{show}}, offering a more detailed
  display of object content.

}
\usage{

detail(object, ...)

}

\arguments{

  \item{object}{An object derived from class \code{ShortRead}. See help
    pages for individual objects, e.g.,
    \code{\linkS4class{ShortReadQ}}. The default simply invokes
    \code{\link{show}}.}

  \item{...}{Additional arguments. The default definition makes no use
    of these arguments.}

}
\value{

  None; the function is invoked for its side effect (detailed display of
  object content).

}
\author{Martin Morgan}
\examples{
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
show(sp)
detail(sp)
}
\keyword{manip}
