\name{RochePath-class}
\docType{class}
\alias{RochePath-class}
\alias{detail,RochePath-method}
\alias{read454,RochePath-method}
\alias{readFastaQual,RochePath-method}
\alias{readFastaQual,character-method}
\alias{readFasta}
\alias{readFasta,character-method}
\alias{readFasta,RochePath-method}
\alias{readPath}
\alias{readQual}
\alias{readQual,character-method}
\alias{readQual,RochePath-method}
\alias{readBaseQuality,RochePath-method}
\alias{read454}
\alias{readFastaQual}
\alias{runNames}
\alias{qualPath}
\alias{runNames,RochePath-method}
\alias{RocheSet,character-method}
\alias{RocheSet,RochePath-method}
\alias{runNames,RochePath-method}
\alias{show,RochePath-method}
%
\alias{RochePath}

\title{"RochePath" class representing a Roche (454) experiment location}

\description{

  This class represents the directory location where Roche (454) result
  files (fasta sequences and qualities) can be found.

}

\section{Objects from the Class}{

  Objects from the class are created with the \code{RochePath}
  constructor:

  \code{RochePath(experimentPath = NA_character_,
    readPath = experimentPath,
    qualPath = readPath, ..., verbose = FALSE)
  }

  \describe{

    \item{experimentPath}{\code{character(1)} or
      \code{\linkS4class{RochePath}} pointing to the top-level directory
      of a Roche experiment.}

    \item{readPath}{\code{character()} of directories (typically in
      \code{experimentPath}) containing sequence (read) information. The
      default selects all directories matching
      \code{list.files(experimentPath, "run")}.}

    \item{qualPath}{\code{character()} of directories (typically in
      \code{experimentPath}) containing quality information. The default
      selects all directories matching \code{list.files(experimentPath,
      "run")}.}

    \item{verbose}{\code{logical(1)} indicating whether invalid paths
      should be reported interactively.}

  }
}
\section{Slots}{
  \code{RocheSet} has the following slots:
  \describe{

    \item{\code{readPath}:}{Object of class \code{"character"}, as
      described in the constructor, above.}

    \item{\code{qualPath}:}{Object of class \code{"character"}, as
      described in the constructor, above.}

    \item{\code{basePath}:}{Object of class \code{"character"},
      containing the \code{experimentPath}.}

  }
}
\section{Extends}{
Class \code{"\linkS4class{ExperimentPath}"}, directly.
Class \code{"\linkS4class{.Roche}"}, directly.
Class \code{"\linkS4class{.ShortReadBase}"}, by class "ExperimentPath", distance 2.
Class \code{"\linkS4class{.ShortReadBase}"}, by class ".Roche", distance 2.
}
\section{Methods}{
  \code{RochePath} has the following methods or functions defined:
  \describe{
    \item{readFasta}{\code{signature(dirPath = "character",
        pattern=".\\.fna$", sample = 1,  ...)}:} 
    \item{readFasta}{\code{signature(dirPath = "RochePath",
        pattern=".\\.fna$", sample = 1, run = 1, ...)}: 

      Read sequences from files matching \code{list.files(dirPath,
      pattern)} (when \code{dirPath="character"}) or
      \code{list.files(readPath(dir)[run], pattern)}, retaining reads
      corresponding to \code{sample}. The result is a
      \code{DNAStringSet}.}

    \item{readQual}{\code{signature(dirPath = "RochePath", reads=NULL, 
        pattern="\\.qual$", sample=1, run=1, ...)}: 

      Read quality scores from files matching
      \code{list.files(qualPath(dirPath)[run])}, corresponding to
      \code{sample}. Non-null \code{reads} is used as an (optional)
      template for parsing quality scores. }

    \item{readFastaQual}{\code{signature(dirPath = "RochePath",
      fastaPattern = "\\.fna$", qualPattern = "\\.qual$", sample = 1,
      run = 1)}: read
      sequences and quality scores into a
      \code{\linkS4class{ShortReadQ}}.}
    \item{readFastaQual}{\code{signature(dirPath = "character",
      fastaPattern = "\\.fna$", qualPattern = "\\.qual$", sample = 1,
      run = 1)}: wrapper for method above, coercing \code{dirPath}
      to a \code{RochePath} via \code{RochePath(dirPath)}.}

    \item{readBaseQuality}{\code{signature(dirPath = "RochePath", ...)}:
      Reads in base and quality information. Currently delegates to
      \code{readFastaQual}, above, but will do more after
      \code{RochePath} supports more file types.
    }
    
    \item{read454}{
      \code{signature(dirPath = "RochePath", ...)}:
      Currently just a wrapper that passes its arguments on to a call to
      \code{readFastaQual}, documented above.
    }
    
    \item{readPath}{\code{signature(object = "RochePath")}: return the
      contents of the \code{readPath} slot.}

    \item{runNames}{\code{signature(object = "RochePath")}: return the
      \code{basename}s of \code{readPath(object)}.}

    \item{RocheSet}{\code{signature(path = "RochePath")}: create a
      \code{\linkS4class{RocheSet}} from \code{path}.}

  }

  Additional methods include:
  \describe{

    \item{show}{\code{signature(object = "RochePath")}: Briefly
      summarize the experiment path locations.}

    \item{detail}{\code{signature(object = "RochePath")}: Provide
      additional detail on the Roche path. All file paths are presented
      in full.}
  }
}

\author{Michael Lawrence <mflawrence@fhcrc.org>}

\seealso{

  \code{\linkS4class{ExperimentPath}}.

}
\examples{
showClass("RochePath")
}
\keyword{classes}
