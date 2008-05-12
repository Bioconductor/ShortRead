\name{SolexaSet-class}
\docType{class}
\alias{SolexaSet-class}
% constructors
\alias{SolexaSet}
\alias{SolexaSet,character-method}
% methods
\alias{detail,SolexaSet-method}
\alias{laneNames,SolexaSet-method}
\alias{laneNames,AnnotatedDataFrame-method}
\alias{show,SolexaSet-method}
% transforming methods
\alias{readAligned,SolexaSet-method}

\title{"SolexaSet" coordinating Solexa output locations with sample annotations}
\description{

  This class coordinates the file hierarchy produced by the Solexa
  `pipeline' with annotation data contained in an
  \code{\linkS4class{AnnotatedDataFrame}} (defined in the \pkg{Biobase}
  package).

}
\section{Objects from the Class}{

  Objects can be created from the constructor:

  \code{SolexaSet(path, ...)}.

  \describe{
    \item{path}{A \code{character(1)} vector giving the fully-qualified
      path to the root of the directory hierarchy associated with each
      Solexa flow cell, or an object of class \code{SolexaPath} (see
      \code{\linkS4class{SolexaPath}} for this method).}

    \item{...}{Additional arguments, especially \code{laneDescription},
	  an \code{\linkS4class{AnnotatedDataFrame}} describing the 
	  content of each of the 8 lanes in the Solexa flow cell.}
  }
}
\section{Slots}{
  \code{SolexaSet} has the following slots:
  \describe{
    \item{\code{solexaPath}:}{Object of class \code{"SolexaPath"}.}
    \item{\code{laneDescription}:}{Object of class
	  \code{"AnnotatedDataFrame"}, containing information about the
	  samples in each lane of the flow cell.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{.Solexa}"}, directly.
Class \code{"\linkS4class{.ShortReadBase}"}, by class ".Solexa", distance 2.
}
\section{Methods}{
  \describe{

    \item{solexaPath}{\code{signature(object = "SolexaSet")}: Return the
      directory paths present when this object was created as a
      \code{\linkS4class{SolexaPath}}.}

    \item{laneNames}{\code{signature(object = "SolexaSet")}: Return the
      names of each lane in the flow cell, currently names are simply
      \code{1:8}. }

    \item{show}{\code{signature(object = "SolexaSet")}: Briefly
      summarize the experiment path and lane description of the Solexa
      set. }

    \item{detail}{\code{signature(object = "SolexaSet")}: Provide
      additional detail on the Solexa set, including the content of
      \code{solexaPath} and the \code{pData} and \code{varMetadata} of
      \code{laneDescription}.}

  }

  Methods transforming \code{SolexaSet} objects include:
  \describe{
    \item{readAligned}{
      \code{signature(dirPath = "SolexaSet", pattern = "s_[1-8]_export.txt", run=1, ...)}: 
      
      Use \code{analysisPath(solexaPath(dirPath))[[run]]} as the
      directory path and \code{pattern="s_[1-8]_export.txt"} as the
      pattern for discovering Eland-aligned reads in the Solexa 'export'
      file format. Note that the default method reads \emph{all} aligned
      read files into a single object; often one will want to specify a
      pattern for each lane.}

  }
}
\author{Martin Morgan}
\examples{
showClass("SolexaSet")
showMethods(class="SolexaSet")
## construct a SolexaSet
sf <- system.file("extdata", package="ShortRead")
df <- data.frame(Sample=c("Sample 1", "Sample 2", "Sample 3", "Sample
                          4", "Center-wide control", "Sample 6", "Sample
                          7", "Sample 8"),

                 Genome=c(rep("hg18", 4), "phi_plus_SNPs.txt",
                          rep("hg18", 3)))
dfMeta <- data.frame(labelDescription=c("Type of sample",
                                        "Alignment genome"))
adf <- new("AnnotatedDataFrame", data=df, varMetadata=dfMeta)
SolexaSet(sf, adf)
}
\keyword{classes}
