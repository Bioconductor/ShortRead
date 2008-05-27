\name{SolexaPath-class}
\docType{class}
\alias{SolexaPath-class}
% constructors
\alias{SolexaPath}
% methods
\alias{SolexaSet,SolexaPath-method}
\alias{show,SolexaPath-method}
\alias{detail,SolexaPath-method}
% transforming methods
\alias{readFastq,SolexaPath-method}
\alias{readAligned,SolexaPath-method}

\title{"SolexaPath" class representing a standard output file hierarchy}
\description{

  Solexa produces a hierarchy of output files. The content of the
  hierarchy varies depending on analysis options. This class represents
  a standard class hierarchy, constructed by searching a 
  file hierarchy for appropriately named directories.

}
\section{Objects from the Class}{

  Objects from the class are created by calls to the constructor:

  \code{SolexaPath(experimentPath,
                   dataPath=.solexaPath(experimentPath, "Data"),
                   scanPath=.solexaPath(dataPath, "GoldCrest"),
                   imageAnalysisPath=.solexaPath(dataPath, "^C"),
                   baseCallPath=.solexaPath(imageAnalysisPath,
                     "^Bustard"),
                   analysisPath=.solexaPath(baseCallPath,
                     "^GERALD"),
                     ..., verbose=FALSE)}

  \describe{
    \item{experimentPath}{\code{character(1)} object pointing to the
      top-level directory of a Solexa run, e.g.,
      \code{/home/solexa/user/080220_HWI-EAS88_0004}. This is the only
      required argument}
    \item{dataPath}{(optional) Solexa \sQuote{Data} folder .}
    \item{scanPath}{(optional) Solexa GoldCrest image scan path.}
    \item{imageAnalysisPath}{(optional) Firecrest image analysis path.}
    \item{baseCallPath}{(optional) Bustard base call path.}
    \item{analysisPath}{(optional) Gerald analysis pipeline path.}
    \item{...}{Additional arguments, unused by currently implemented
      methods.}
    \item{verbose=FALSE}{(optional) logical vector which, when
      \code{TRUE} results in warnings if paths do not exist.}
  }
  All paths must be fully-specified.
}
\section{Slots}{
  \code{SolexaPath} has the following slots, containing either a fully
  specified path to the corresponding directory (described above) or
  \code{NA} if no appropriate directory was discovered.
  \describe{
    \item{\code{experimentPath}}{See above.}
    \item{\code{dataPath}}{See above.}
    \item{\code{scanPath}}{See above.}
    \item{\code{imageAnalysisPath}}{See above.}
    \item{\code{baseCallPath}}{See above.}
    \item{\code{analysisPath}}{See above.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{.Solexa}"}, directly.
Class \code{"\linkS4class{.ShortReadBase}"}, by class ".Solexa", distance 2.
}
\section{Methods}{
  Transforming methods include:
  \describe{

    \item{readFastq}{

      \code{signature(dirPath = "SolexaPath", pattern = "s_[1-8]_sequence.txt", run = 1, ...)}:

      Use \code{analysisPath(dirPath)[[run]]} as the directory path and
      \code{pattern="s_[1-8]_sequence.txt"} as the pattern for
      discovering fastq-formatted files, returning a
      \code{\linkS4class{ShortReadQ}} object. Note that the default
      method reads \emph{all} sequence files into a single object; often
      one will want to specify a pattern for each lane. }

    \item{readAligned}{
      \code{signature(dirPath = "SolexaPath", pattern = "s_[1-8]_export.txt", run=1, ...)}: 
      
      Use \code{analysisPath(dirPath)[[run]]} as the directory path and
      \code{pattern="s_[1-8]_export.txt"} as the pattern for discovering
      Eland-aligned reads in the Solexa 'export' file format. Note that
      the default method reads \emph{all} aligned read files into a
      single object; often one will want to specify a pattern for each
      lane.}

    \item{SolexaSet}{\code{signature(path = "SolexaPath")}: create a
      \code{\linkS4class{SolexaSet}} object based on \code{path}.}

  }

  Additional methods include:
  \describe{

    \item{show}{\code{signature(object = "SolexaPath")}: briefly
      summarize the file paths of \code{object}. The
      \code{experimentPath} is given in full; the remaining paths are
      identified by their leading characters.}

    \item{detail}{\code{signature(object = "SolexaPath")}: summarize
      file paths of \code{object}. All file paths are presented in
      full.}
  }
}

\author{Martin Morgan}
\examples{
showClass("SolexaPath")
showMethods(class="SolexaPath")
sf <- system.file("extdata", package="ShortRead")
sp <- SolexaPath(sf)
sp
readFastq(sp, pattern="s_1_sequence.txt")
}
\keyword{classes}
