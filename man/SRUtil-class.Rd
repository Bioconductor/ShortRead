\name{SRUtil-class}
\docType{class}

\alias{.SRUtil-class}
\alias{SRError-class}
\alias{SRWarn-class}
\alias{SRList-class}
\alias{SRVector-class}
% constructors
\alias{SRError}
\alias{SRWarn}
\alias{SRList}
\alias{SRVector}
% methods
\alias{detail,.ShortReadBase-method} % internal method, not documented
\alias{length,SRList-method}
\alias{names,SRList-method}
\alias{names<-,SRList,character-method}
\alias{lapply,SRList-method}
\alias{lapply,SRList,ANY-method}
\alias{sapply,SRList-method}
\alias{[,SRList,ANY,missing-method}
\alias{[,SRList,ANY,missing,ANY-method}
\alias{[[,SRList,ANY,missing-method}
\alias{detail,SRList-method}
\alias{show,SRList-method}

\alias{detail,SRVector-method}
\alias{show,SRVector-method}


\title{".SRUtil" and related classes}
\description{

  These classes provide important utility functions in the
  \pkg{ShortRead} package, but may occasionally be seen by the user and
  are documented here for that reason.

}
\section{Objects from the Class}{

  Utility classes include:

  \itemize{

    \item \code{.SRUtil-class} a virtual base class from which all
    utility classes are derived.

    \item \code{SRError-class} created when errors occur in
    \pkg{ShortRead} package code.

    \item \code{SRWarn-class} created when warnings occur in
    \pkg{ShortRead} package code

    \item \code{SRList-class} representing a list (heterogeneous
    collection) of objects.

    \item \code{SRVector-class} representing a vector (homogeneous
    collection, i.e., all elements of the same class) of objects.

  }

  Objects from these classes are not normally constructed by the
  user. However, constructors are available, as follows.

  \code{SRError(type, fmt, ...)}, \code{SRWarn(type, fmt, ...)}:
  \describe{

    \item{type}{\code{character(1)} vector describing the type of the
      error. \code{type} must come from a pre-defined list of types.}

    \item{fmt}{a \code{\link{sprintf}}-style format string for the
      message to be reported with the error.}

    \item{...}{additional arguments to be interpolated into \code{fmt}.}
  
  }

  \code{SRList(...)}
  \describe{

    \item{...}{elements of any type or length to be placed into the
      \code{SRList}. If the length of \code{...} is 1 and the 
      argument is a list, then the list itself is placed into
      \code{SRList}.}

  }

  \code{SRVector(..., vclass)}
  \describe{

    \item{...}{elements all satisfying an \code{\link{is}} relationship
      with \code{vclass}, to be placed in \code{SRVector}.}

    \item{vclass}{the class to which all elements in \code{...}
      belong. If \code{vclass} is missing and \code{length(list(...))}
      is greater than zero, then \code{vclass} is taken to be the class
      of the first argument of \code{...}.}

  }
  \code{SRVector} errors:
  \describe{

    \item{SRVectorClassDisagreement}{this error occurs when not all
    arguments \code{...} satisfy an \sQuote{is} relationship with
    \code{vclass}.}

  }

}
\section{Slots}{
  \code{SRError} and \code{SRWarn} have the following slots defined:
  \describe{

    \item{\code{.type}:}{Object of class \code{"character"} containing
      the type of error or warning. \code{.type} must come from a
      pre-defined list of types, see, e.g.,
      \code{ShortRead:::.SRError_types}.}

    \item{\code{.message}:}{Object of class \code{"character"}
      containing a detailed message describing the error or warning.}

  }

  \code{SRList} has the following slot defined:
  \describe{

    \item{\code{.srlist}:}{Object of class \code{"list"} containing the
      elements in the list.}
    
  }

  \code{SRVector} extends \code{SRList}, with the following additional
  slot:
  \describe{

    \item{\code{vclass}:}{Object of class \code{"character"} naming the
      type of object all elements of \code{SRVector} must be.}

  }
}
\section{Methods}{

  Accessors are available for all slots, and have the same name as the
  slot, e.g., \code{vclass} to access the \code{vclass} slot of
  \code{SRVector}.  Internal slots (those starting with \sQuote{.} also
  have accessors, but these are not exported e.g.,
  \code{ShortRead:::.type}.

  \code{SRList} has the following methods:

  \describe{

    \item{length}{\code{signature(x = "SRList")}: return the
      (\code{integer(1)}) length of the \code{SRList}.}

    \item{names}{\code{signature(x = "SRList")}: return a character
      vector of list element names. The length of the returned vector is
      the same as the length of \code{x}.}

    \item{names<-}{\code{signature(x = "SRList", value = "character")}:
      assign \code{value} as names for members of \code{x}.}

    \item{[}{\code{signature(x = "SRList", i = "ANY", j = "missing")}:
      subset the list using standard R list subset paradigms.}

    \item{[[}{\code{signature(x = "SRList", i = "ANY", j = "missing")}:
      select element \sQuote{i} from the list, using standard R list
      selection paradigms.}

    \item{lapply}{\code{signature(X = "SRList", FUN="ANY")}: apply a function to
      all elements of \code{X}, with additional arguments interpreted
      as with \code{\link{lapply}}.}

    \item{sapply}{\code{signature(X = "SRList")}: apply a function to
      all elements of \code{X}, simplifying the result if
      possible. Additional arguments interpreted as with
      \code{\link{sapply}}.}

    \item{show}{\code{signature(object = "SRList")}: display an
      informative summary of the object content, including the length of
      the list represented by \code{object}.}

    \item{detail}{\code{signature(x = "SRList")}: display a more
      extensive version of the object, as one might expect from printing
      a standard list in R.}

  }

  \code{SRVector} inherits all methods from \code{SRList}, and has the
  following additional methods:

  \describe{

    \item{show}{\code{signature(object = "SRVector")}: display an
      informative summary of the object content, e.g., the vector class
      (\code{vclass}) and length.}

    \item{detail}{\code{signature(x = "SRVector")}: display a more
      extensive version of the object, as one might expect from a
      printing a standard R list.}

  }
}

\author{Martin Morgan}
\examples{
getClass(".SRUtil", where=getNamespace("ShortRead"))
ShortRead:::.SRError_types
ShortRead:::.SRWarn_types

detail(SRList(1:5, letters[1:5]))

tryCatch(SRVector(1:5, letters[1:5]),
         SRVectorClassDisagreement=function(err) {
             cat("caught:", conditionMessage(err), "\n")
         })
}
\keyword{classes}
