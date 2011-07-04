\name{Snapshot-class}

\Rdversion{1.1}
\docType{class}

% Constructor:
\alias{Snapshot}
\alias{Snapshot,character,GRanges-method}
\alias{Snapshot,character,missing-method}

% Accessors:
\alias{files}
\alias{files,Snapshot-method}
\alias{functions}
\alias{functions,Snapshot-method}
\alias{show,Snapshot-method}
\alias{view}
\alias{view,Snapshot-method}
\alias{vrange}
\alias{vrange,Snapshot-method}

% methods: 
\alias{pan}
\alias{pan,Snapshot-method}
\alias{togglefun}
\alias{togglefun,Snapshot-method}
\alias{togglep}
\alias{togglep,Snapshot-method}
\alias{togglez}
\alias{togglez,Snapshot-method}
\alias{zoom}
\alias{zoom,Snapshot-method}

\title{Class \code{"Snapshot"}}

\description{

  A \code{\linkS4class{Snapshot}}-class to visualize genomic data from 
  BAM files with zoom and pan functionalities.

}

\usage{
Snapshot(files, range, ...)
}

\arguments{

  \item{files}{A character() specifying the file(s) to be visualized.}

  \item{range}{A \code{\link{GRanges}} object specifing the range to be
    visualized.}

   \item{...}{Additional, optional, arguments to be passed to the
     Snapshot \code{initialize} function. Arguments include:

     \describe{

       \item{functions:}{A \code{\link{SnapshotFunctionList}} of
	 functions, in addition to built-in \sQuote{fine_coverage},
	 \sQuote{coarse_coverage}, \sQuote{multifine_coverage}, to be
	 used for visualization.}

       \item{currentFunction:}{character(1) naming the function, from
	 \code{functinos} to be used for data input and
	 visualization. The default chooses a function based on the
	 scale at which the data is being visualized.}

       \item{.auto_display:}{logical(1) indicating whether the
	 visualization is to be updated when \code{show} is invoked.}

       \item{.debug}{logical(1) indicating whether debug messages are to
	 be printed.}

   }}

}


\section{Constructor}{
  \describe{
    \item{}{
       \code{Snapshot(files, param, ..., functions,
       currentFunction)}: Creates a Snapshot object. If
       \code{currentFunction} is not specified, default \code{functions}
       are used to evaluate imported recored and to view the coverage
       on selected range.

       \describe{
          \item{\code{files}:}{A character vector naming the BAM files
             to be imported. If multiple files, sum of coverage is
             returned (see \code{.data} field)  and visualized.}
          \item{\code{params}:}{A \code{GRanges} or missing object
             specifying which records are to be imported and viewed.}
          \item{\code{functions}:}{A \code{SnapshotFunction} object
             which records are to be evaluated and visualized.}
          \item{\code{currentFunction}:}{A character vector naming the
             function that the records are to be immediately evaluated.}
          \item{\code{\ldots}:}{Optional arguments.}
       }
    }
  }
}

\section{Methods}{
  \describe{
    \item{zoom}{\code{signature(x = "Snapshot")}: Zoom (in or out) the
       current plot. } 
    \item{pan}{\code{signature(x = "Snapshot")}: Pan (right or left)
       the current plot. } 
    \item{togglefun}{\code{signature(x = "Snapshot")}: Toggle the current
       functions which imported records are to be immediately evaluated.  }
    \item{togglep}{\code{signature(x = "Snapshot")}: Toggle the panning
       effects.}
    \item{togglez}{\code{signature(x = "Snapshot")}: Toggle the zooming 
       effects.} 
  }
}

\section{Accessors}{
   \describe{
     \item{show}{\code{signature(object = "Snapshot")}: Display a
       \code{Snapshot} object. }

     \item{files}{\code{signature(x = "Snapshot")}: Get the \code{files}
       field (object of class \code{BamFileList}) of a \code{Snapshot} object.}

     \item{functions}{\code{signature(x = "Snapshot")}: Get the
       \code{functions} field (object of \code{SnapshotFunctionList}) of a
       \code{Snapshot} object.}

     \item{view}{\code{signature(x = "Snapshot")}: Get the \code{view}
       field (object of \code{SpTrellis}) of a \code{Snapshot} object.}

     \item{vrange}{\code{signature(x = "Snapshot")}: Get the
       \code{.range} field (object of \code{GRanges}) of a \code{Snapshot}
       object. }
  }
}

\section{Fields}{
  \describe{
    \item{\code{.debug}:}{Object of class \code{function} to display
    messages while in debug mode }
    \item{\code{.auto_display}:}{Object of class \code{logical} to
       automatically disply the coverage plot. }
    \item{\code{.range}:}{Object of class \code{GRanges} indicating
      which ranges of records to be imported from BAM fiels. }
    \item{\code{.zin}:}{Object of class \code{logical} indicating
      whether the current zooming effect is zoom in. }
    \item{\code{.pright}:}{Object of class \code{logical} indicating
      whether the current panning effect is right. }
    \item{\code{.data}:}{Object of class \code{data.frame} containing
      coverage a position is represented for each strand and BAM file.}
    \item{\code{.data_dirty}:}{Object of class \code{logical}
      indicating whether to re-evaluate the imported records. }
    \item{\code{.initial_functions}:}{Object of class 
       \code{SnapshotFunctionList} available by the \code{Snapshot} object. }
    \item{\code{.current_function}:}{Object of class \code{character}
      of the function the imported recored are currently evaluated and
      visualized.}
    \item{\code{functions}:}{Object of class
      \code{SnapshotFunctionList} of customized functions to evaluate
       and visualize the imported records. }
    \item{\code{files}:}{Object of class \code{BamFileList} to be imported. }
    \item{\code{view}:}{Object of class \code{SpTrellis} that is
    essentially a reference class wrapper of \code{Trellis} objects.  }
  }
}

\section{Class-Based Methods}{
  \describe{
    \item{\code{display()}:}{Display the current \code{Snapshot} object. }
    \item{\code{pan()}:}{Pan (right or left) the current plot. }
    \item{\code{zoom()}:}{Zoom (in or out) the current plot. }
    % \item{\code{set_range(range)}:}{ ~~ }
    \item{\code{toggle(zoom, pan, currentFunction)}:}{Toggle zooming,
    panning effects or the currentFuction in which the imported
    records are to be evaluated and visualized.}
    % \item{\code{initialize(..., functions, currentFunction, .range, .auto_display, .debug)}:}{ ~~ }
  }
}


\author{Chao-Jen \url{cwon2@fhcrc.org}}

\seealso{\code{\link{SpTrellis}}}

\examples{
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
s <- Snapshot(fl)

}

\keyword{classes}
