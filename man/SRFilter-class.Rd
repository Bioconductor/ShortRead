\name{SRFilter-class}
\docType{class}
\alias{SRFilter-class}
\alias{srFilter,SRFilter-method}
\alias{name}
\alias{name,SRFilter-method}
\alias{show,SRFilter-method}

\title{"SRFilter" for representing functions operating on ShortRead objects}
\description{

  Objects of this class are functions that, when provided an appropriate
  object from the ShortRead package, return logical vectors indicating
  which parts of the object satisfy the filter criterion.

  A number of filters are built-in (described below); users are free to
  create their own filters, using the \code{srFilter} function.

}

\section{Objects from the Class}{

  Objects can be created through \code{\link{srFilter}} (to create a
  user-defined filter) or through calls to constructors for predefined
  filters, as described on the \code{\link{srFilter}} page.

}

\section{Slots}{
  \describe{

    \item{\code{.Data}:}{Object of class \code{"function"} taking a
      single named argument \code{x} corresponding to the ShortRead
      object that the filter will be applied to. The return value of the
      filter function is expected to be a logical vector that can be
      used to subset \code{x} to include those elements of \code{x}
      satisfying the filter.
    }

    \item{\code{name}:}{Object of class \code{"ScalarCharacter"}
      representing the name of the filter. The name is useful for
      suggesting the purpose of the filter, and for debugging failed
      filters.
    }
  }
}

\section{Extends}{
Class \code{"\linkS4class{function}"}, from data part.
Class \code{"\linkS4class{.SRUtil}"}, directly.
Class \code{"\linkS4class{OptionalFunction}"}, by class "function", distance 2.
Class \code{"\linkS4class{PossibleMethod}"}, by class "function", distance 2.
}
\section{Methods}{
  \describe{

    \item{srFilter}{\code{signature(fun = "SRFilter")}: Return the function
      representing the underlying filter; this is primarily for
      interactive use to understanding filter function; usually the
      filter is invoked as a normal function call, as illustrated below
    }

    \item{name}{\code{signature(x = "SRFilter")}: Return, as a
      \code{ScalarCharacter}, the name of the function.
    }

    \item{show}{\code{signature(object = "SRFilter")}: display a brief
      summary of the filter
    }

  }
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{
  \code{\link{srFilter}} for predefined and user-defined filters.
}

\examples{
## see ?srFilter
}

\keyword{classes}
