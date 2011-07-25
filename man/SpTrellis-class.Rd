\name{SpTrellis-class}

\Rdversion{1.1}
\docType{class}

% constructor:
\alias{SpTrellis-class}
\alias{SpTrellis}

% methods:
\alias{zi, SpTrellis-method}
\alias{zo, SpTrellis-method}
\alias{right, SpTrellis-method}
\alias{left, SpTrellis-method}
\alias{restore, SpTrellis-method}
\alias{show, Sptrellis-method}

\title{Class "SpTrellis"}

\description{

  A reference class to manage the trellis graphics related component of the 
  \code{\link{Snapshot}} functionality for visualization of genomic data.
 
}

\usage{
SpTrellis(trellis, debug_enabled=FALSE)
}

\arguments{
  \item{trellis}{A trellis object for storing the plot of the genome area being 
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

\section{Methods}{
  \describe{
    \item{zi}{\code{signature(x="SpTrellis")}: zoom in}
    \item{zo}{\code{signature(x="SpTrellis")}: zoom out}
    \item{right}{\code{signature(x="SpTrellis")}: shift to the right}
    \item{left}{\code{signature(x="SpTrellis")}: shift to the left}
    \item{restore}{\code{signature(x="SpTrellis")}: restore to the
      original plot}
    \item{show}{\code{signature(x="SpTrellis")}: show the current plot}
    \item{update}{\code{signature(x="SpTrellis")}: update the trellis
      parameters of the \code{SpTrellis} object.}
    
  }
}

\author{Chao-Jen \url{cwon2@fhcrc.org}}

\seealso{\code{\link{Snapshot}}}

\examples{
data(yeastChrI) # data frame
col <- c("#66C2A5", "#FC8D62")
cv <- xyplot(data ~ pos, data=yeastChrI, group=group, type="s",
             col=col, main="yeast chrI:1 - 2e5",
             ylab="Coverage", xlab="Coordinate",
             scales=list(y=list(tck=c(1,0)),
                             x=list(rot=45, tck=c(1,0), tick.number=20)),
             panel=function(...) {
                     panel.xyplot(...)
                     panel.grid(h=-1, v=20)
                     panel.abline(a=0, b=0, col="grey")
             })
s <- SpTrellis(cv)
s
zi(s)
zi(s)
right(s)
left(s)
zo(s)
restore(s)

}

\keyword{classes}
