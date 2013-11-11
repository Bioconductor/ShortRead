\name{srapply}
\alias{srapply}
\title{DEPRECATED Apply-like function for distribution across MPI-based clusters.}
\description{

  DEPRECATED; use the BiocParallel package instead.

  This \code{lapply} like function evaluates locally or, if \pkg{Rmpi}
  or \pkg{parallel} is loaded (and \pkg{Rmpi} workers spawned) and
  \code{options(srapply_fapply)} set appropriately (see below), across
  nodes in a cluster. Errors in evaluation of \code{FUN} generate
  warnings; results are trimmed to exclude results where the error
  occurs.

}

\usage{
srapply(X, FUN, ..., fapply = .fapply(), reduce = .reduce(),
        USE.NAMES = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{X}{Tasks to be distributed. \code{X} should be an object for
    which \code{lapply} or \code{sapply} are defined (more precisely,
    \code{mpi.parLapply}, \code{mpi.parSapply}, or
    \code{mclapply}). Performance is best when these objects are
    relatively small, e.g., file names, compared to the work to be done
    on each by \code{FUN}.}

  \item{FUN}{A function to be applied to each element of \code{X}. The
    function must have \code{...} or named argument \code{verbose} in
    its signature. It is best if it makes no reference to variables
    other than those in its argument list. or in loaded packages (the
    \pkg{ShortRead} package is available on remote nodes).}

  \item{\dots}{Additional arguments, passed to \code{FUN}.}

  \item{fapply}{An optional argument defining an \code{lapply}-like
    function to be used in partitioning \code{X}. See details, below.}

  \item{reduce}{Optional function accepting a list (the result of
    \code{fapply} and summarizing this. The default reports errors in
    function evaluation as warnings, returning the remaining values as
    elements of a list. See details below for additional hints.}

  \item{USE.NAMES}{If \code{TRUE} and if \code{X} is character,
    use \code{X} as \code{names} for the result unless it had names
    already.}

  \item{verbose}{Report whether evaluation is local or mpi-based; also
    forwarded to \code{FUN}, allowing detailed reports from remote
    instances.}

}

\details{

    The default value for \code{fapply} is available with
    \code{ShortRead:::.fapply()}. It tests
    \code{getOption("srapply_fapply")} for value \dQuote{Rmpi} or
    \dQuote{parallel}.

    If \pkg{Rmpi} is indicated, \code{fapply} ensures that
    \pkg{ShortRead} is \code{require}d on all workers, and then invokes
    \code{mpi.parLapply} with arguments \code{X}, \code{FUN},
    \code{...}, and \code{verbose}. The function \code{FUN} is wrapped
    so that errors are returned as objects of class \code{SRError} with
    type \code{RemoteError}.  If no workers are available, the code
    evaluates \code{FUN} so that errors are reported as with remote
    evaluation.

    If \pkg{parallel} is indicated, \code{fapply} invokes
    \code{mclapply} with arguments as for \code{mpi.parLapply}.

    Custom \code{reduce} functions might be written as
    \code{reduce=function(lst) unlist(lst, use.names=TRUE)}.

}
\value{

  The returned value depends on the value of \code{reduce}, but by
  default is a list with elements containing the results of \code{FUN}
  applied to each of \code{X}. Evaluations resulting in an error have
  been removed, and a warning generated.

}
\author{Martin Morgan <mtmorgan@fhcrc.org>}

\examples{
\dontrun{
## DEPRECATED; use BiocParallel::bplapply instead
srapply(1:10, function(i, ...) i)
## collapse result to vector
srapply(1:10, function(i, ...) i, reduce=unlist)
x <- srapply(1:10, function(i, ...) {
    if (runif(1)<.2) stop("oops") else i
})
length(x) ## trimmed to exclude errors
}
}
\keyword{manip}
