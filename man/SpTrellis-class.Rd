\name{SpTrellis-class}

\Rdversion{1.1}
\docType{class}

\alias{SpTrellis-class}
\alias{SpTrellis}

\title{Class "SpTrellis"}

\description{

  A reference class to manage the trellis graphics related component of the 
  \code{\link{Snapshot}} functionality for visualization of genomic data.
 
}

\usage{
SpTrellis(trellis, debug_enabled=FALSE)
}

\arguments{
  \item{trellis}{ A trellis device for storing the plot of the genome area being 
      visualized.}
  \item{debug_enabled}{\code{logical(1)} indicating whether class methods
          should report debugging information to the user.}
}

\section{Fields}{
  \describe{
    \item{\code{trellis}:}{Object of class \code{trellis} for storing the plot
        information.}
    \item{debug_enabled}{\code{logical(1)} indicating whether class methods
          should report debugging information to the user.}
  }
}

\author{Chao-Jen \url{cwon2@fhcrc.org}}

\seealso{\code{\link{Snapshot}}}

\examples{showClass("SpTrellis")}

\keyword{classes}
