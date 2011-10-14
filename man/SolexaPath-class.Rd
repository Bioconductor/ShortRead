\name{SolexaPath-class}
\docType{class}
\alias{SolexaPath-class}
% constructors
\alias{SolexaPath}
% methods
\alias{SolexaSet,SolexaPath-method}
\alias{qa,SolexaPath-method}
\alias{report,SolexaPath-method}
\alias{show,SolexaPath-method}
\alias{detail,SolexaPath-method}
% transforming methods
\alias{readIntensities,SolexaPath-method}
\alias{readPrb,SolexaPath-method}
\alias{readQseq,SolexaPath-method}
\alias{readFasta,SolexaPath-method}
\alias{readFastq,SolexaPath-method}
\alias{readBaseQuality,SolexaPath-method}
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
                   imageAnalysisPath=.solexaPath(dataPath, "^(C|IPAR)"),
                   baseCallPath=.solexaPath(imageAnalysisPath,
                     "^Bustard"),
                   analysisPath=.solexaPath(baseCallPath,
                     "^GERALD"),
                   ..., verbose=FALSE)
		 }

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
    \item{\code{basePath}}{See \code{experimentPath}, above.}
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

    \item{readIntensities}{

      \code{signature(dirPath = "SolexaPath", pattern=character(0), run, ...)}:

      Use \code{imageAnalysisPath(sp)[run]} as the directory path(s) and
      \code{pattern=character(0)} as the pattern for discovering Solexa
      intensity files. See
      \code{\link{readIntensities,character-method}} for additional
      parameters.}
    
    \item{readPrb}{

      \code{signature(dirPath = "SolexaPath", pattern=character(0), run, ...)}:

      Use \code{baseCallPath(dirPath)[run]} as the directory path(s) and
      \code{pattern=character(0)} as the pattern for discovering Solexa
      \sQuote{prb} files, returning a \code{\linkS4class{SFastqQuality}}
      object containing the maximum qualities found for each base of
      each cycle.

      The \code{...} argument may include the named argument
      \code{as}. This influences the return value, as explained on the
      \code{\link{readPrb,character-method}} page.

    }

    \item{readFasta}{

      \code{signature(dirPath, pattern = character(0), ...,
        nrec=-1L, skip=0L)}:

      Use \code{analysisPath(dirPath)[run]} as the directory path(s) for
      discovering fasta-formatted files, returning a
      \code{\linkS4class{ShortRead}} object. The default method reads
      \emph{all} files into a single object.}

    \item{readFastq}{

      \code{signature(dirPath = "SolexaPath", pattern = ".*_sequence.txt",
		run, ..., qualityType="SFastqQuality")}:

      Use \code{analysisPath(dirPath)[run]} as the directory path(s) and
      \code{pattern=".*_sequence.txt"} as the pattern for discovering
      fastq-formatted files, returning a \code{\linkS4class{ShortReadQ}}
      object. The default method reads \emph{all} sequence files into a
      single object.}

    \item{readBaseQuality}{
      \code{signature(dirPath = "SolexaPath", seqPattern = ".*_seq.txt", prbPattern = "s_[1-8]_prb.txt", run, ...)}: 
      
      Use \code{baseCallPath(dirPath)[run]} as the directory path(s) and
      \code{seqPattern=".*_seq.txt"} as the pattern for discovering
      base calls and \code{prbPattern=".*_prb.txt"} as the pattern
      for discovering quality scores. Note that the default method reads
      \emph{all} base call and quality score files into a single object;
      often one will want to specify a pattern for each lane.}

	\item{readQseq}{
	  \code{signature(directory="SolexaPath", pattern=".*_qseq.txt.*", run, ...., filtered=FALSE)}:

      Use \code{analysisPath(dirPath)[run]} as the directory path and
      \code{pattern=".*_qseq.txt.*"} as the pattern for discovering read
      and quality scores in Solexa 'qseq' files. Data from \emph{all}
      files are read into a single object; often one will want to
      specify a pattern for each lane. Details are as for
      \code{\link{readQseq,character-method}}.}

    \item{readAligned}{
      \code{signature(dirPath = "SolexaPath", pattern = ".*_export.txt.*", run, ..., filter=srFilter())}: 
      
      Use \code{analysisPath(dirPath)[run]} as the directory path and
      \code{pattern=".*_export.txt"} as the pattern for discovering
      Eland-aligned reads in the Solexa 'export' file format. Note that
      the default method reads \emph{all} aligned read files into a
      single object; often one will want to specify a pattern for each
      lane. Use an object of \code{\linkS4class{SRFilter}} to select
      specific chromosomes, strands, etc.}

    \item{qa}{
      \code{signature(dirPath="SolexaPath", pattern="character(0)", run, ...)}:

      Use \code{analysisPath(dirPath)[run]} as the directory path(s) and
      \code{pattern=".*_export.txt"} as the pattern for discovering
      Solexa \code{export}-formatted fileds, returning a
      \code{\linkS4class{SolexaExportQA}} object summarizing quality
      assessment. If \code{Rmpi} or \code{parallel} has been initiated,
      quality assessment calculations are distributed across available
      nodes or cores (one node per export file.)}

    \item{report}{
      \code{signature(x, ..., dest=tempfile(), type="pdf")}: Use
      \code{qa(x, ...)} to generate quality assessment measures, and
      use these to generate a quality assessment report at location
      \code{dest} of type \code{type} (e.g., \sQuote{pdf}).
    }

    \item{SolexaSet}{\code{signature(path = "SolexaPath")}: create a
      \code{\linkS4class{SolexaSet}} object based on \code{path}.}

  }

  Additional methods include:
  \describe{

    \item{show}{\code{signature(object = "SolexaPath")}: briefly
      summarize the file paths of \code{object}. The
      \code{experimentPath} is given in full; the remaining paths are
      identified by their leading characters.}

    \item{detail}{\code{signature(x = "SolexaPath")}: summarize
      file paths of \code{x}. All file paths are presented in
      full.}
  }
}

\author{Martin Morgan}
\examples{
showClass("SolexaPath")
showMethods(class="SolexaPath", where=getNamespace("ShortRead"))
sf <- system.file("extdata", package="ShortRead")
sp <- SolexaPath(sf)
sp
readFastq(sp, pattern="s_1_sequence.txt")
\dontrun{
nfiles <- length(list.files(analysisPath(sp), "s_[1-8]_export.txt"))
library(Rmpi)
mpi.spawn.Rslaves(nslaves=nfiles)
report(qa(sp))
}
\dontrun{
nfiles <- length(list.files(analysisPath(sp), "s_[1-8]_export.txt"))
library(parallel)
options(srapply_fapply="parallel")
report(qa(sp))
}
}
\keyword{classes}
