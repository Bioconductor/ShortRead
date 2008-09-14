\name{ExperimentPath-class}
\docType{class}
\alias{ExperimentPath-class}
% constructors
\alias{ExperimentPath}
% etc
\alias{show,ExperimentPath-method}
\alias{detail,ExperimentPath-method}

\title{"ExperimentPath" class representing a file hierarchy of data
  files}

\description{

  Short read technologies often produce a hierarchy of output files. The
  content of the hierarchy varies. This class represents the root of the
  file hierarchy. Specific classes (e.g.,
  \code{\linkS4class{SolexaPath}}) represent different technologies.

}
\section{Objects from the Class}{

  Objects from the class are created by calls to the constructor:
  
  \code{ExperimentPath(experimentPath)}

  \describe{
    \item{experimentPath}{\code{character(1)} object pointing to the
      top-level directory of the experiment; see specific technology
      classes for additional detail.}
    \item{verbose=FALSE}{(optional) logical vector which, when
      \code{TRUE} results in warnings if paths do not exist.}
  }
  All paths must be fully-specified.
}
\section{Slots}{
  \code{ExperimentPath} has one slot, containing a fully
  specified path to the corresponding directory (described above).
  \describe{
    \item{\code{basePath}}{See above.}
  }
  The slot is accessed with \code{experimentPath}.
}
\section{Extends}{
  Class \code{"\linkS4class{.ShortReadBase}"}, directly.
}
\section{Methods}{
  Methods include:
  \describe{

    \item{show}{\code{signature(object = "ExperimentPath")}: briefly
      summarize the file paths of \code{object}.}

    \item{detail}{\code{signature(object = "ExperimentPath")}: summarize
      file paths of \code{object}.}
  }
}

\author{Michael Lawrence}
\examples{
showClass("ExperimentPath")
}
\keyword{classes}
