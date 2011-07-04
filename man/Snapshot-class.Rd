\name{Snapshot-class}

\Rdversion{1.1}
\docType{class}

\alias{Snapshot-class}
\alias{Snapshot}
\alias{Snapshot,character,GRanges-method}
\alias{Snapshot,character,missing-method}
\alias{show,Snapshot-method}
\alias{vrange}
\alias{vrange,Snapshot-method}
\title{Class "Snapshot"}

\description{

  A \code{\linkS4class{Snapshot}}-class to visualize genomic data from 
  BAM files.

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

\section{Fields}{
  Fields are documented for developer use. Fields starting with
  \sQuote{.} are meant to be private (used internally only)

  \describe{
    \item{\code{.debug}:}{Object of class \code{function}.}
    \item{\code{.auto_display}:}{Object of class \code{logical}.}
    \item{\code{.range}:}{Object of class \code{GRanges}.}
    \item{\code{.zin}:}{Object of class \code{logical}.}
    \item{\code{.pright}:}{Object of class \code{logical}.}
    \item{\code{.data}:}{Object of class \code{data.frame}. A data.frame
      holding the Genomic data to be visualized.}
    \item{\code{.data_dirty}:}{Object of class \code{logical}.}
    \item{\code{.initial_functions}:}{Object of class
      \code{SnapshotFunctionList}. A list of functions that are used for
      reading data from the file}
    \item{\code{.current_function}:}{Object of class
      \code{character}. Name of the current function being used by the
      file to read data.}
    \item{\code{functions}:}{Object of class
      \code{SnapshotFunctionList}. A list of functions that are used for
      reading data.}
   \item{\code{files}:}{Object of class \code{BamFileList}. Bam files
     from which data is to be read}
    \item{\code{view}:}{Object of class \code{SpTrellis}. A trellis
      object for storing the current plot information.}
  }
}

\section{Class-Based Methods}{
  \describe{
    \item{\code{initialize(...,functions, currentFunction, .range,
        .auto_display, .debug)}:}{ Called during object creation with values
      to initialize fields}

    \item{\code{set_range(range)}:}{.}
    \item{\code{display()}:}{.}
    \item{\code{.initialize_currentFunction()}:}{.}
    \item{\code{.change_current_function()}:}{.}
    \item{\code{.check_currentFunction()}:}{.}
    \item{\code{.update_data()}:}{.}
    \item{\code{.update_range()}:}{.}
    \item{\code{.initial_range()}:}{.}
    \item{\code{.stop()}:}{.}
    \item{\code{.message()}:}{.}
  }
}

\author{Chao-Jen \url{cwon2@fhcrc.org}}

\seealso{\code{\link{SpTrellis}}}

\examples{
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
s <- Snapshot(fl)

}

\keyword{classes}
