\name{QA-class}
\docType{class}
\alias{QA-class}

\alias{.QA2-class}
\alias{QA-class}
\alias{QAAdapterContamination-class}
\alias{QAFrequentSequence-class}
\alias{QANucleotideByCycle-class}
\alias{QANucleotideUse-class}
\alias{QAQualityByCycle-class}
\alias{QAQualityUse-class}
\alias{QAReadQuality-class}
\alias{QASequenceUse-class}

\alias{QASource-class}
\alias{QAFastqSource-class}
\alias{QAData-class}

\alias{QACollate-class}

\alias{QAFiltered-class}
\alias{QAFlagged-class}
\alias{QASummary-class}

\title{(Updated) classes for representing quality assessment results}
\description{

  Classes derived from \code{.QA-class} represent results of quality
  assurance analyses.

}

\section{Objects from the Class}{

  Users create instances of many of these classes by calling the
  corresponding constructors, as documented on the help page for
  \code{\link{qa2}}. Classes constructed in this way include
  \code{\link{QACollate}}, \code{\link{QAFastqSource}},
  \code{\link{QAAdapterContamination}},
  \code{\link{QAFrequentSequence}}, \code{\link{QANucleotideByCycle}},
  \code{\link{QANucleotideUse}}, \code{\link{QAQualityByCycle}},
  \code{\link{QAQualityUse}}, \code{\link{QAReadQuality}}, and
  \code{\link{QASequenceUse}}.

  The classes \code{QASource}, \code{QAFiltered}, \code{QAFlagged} and
  \code{QASummary} are generated internally, not by users.

}

\section{Extends}{

  \code{.QA2} extends class \code{"\linkS4class{.ShortReadBase}"},
  directly.

  \code{QASummary} is a virtual class extending \code{.QA2}; all
  user-creatable classes extend \code{QASummary}.

  \code{QASource} extends \code{QASummary}. All classes used to
  represent raw data input (\code{QAFastqSource}) extend
  \code{QASource}.

  \code{QAData} is a reference class, used to contain a single instance
  of the fastq used in all QA Summary steps.

  \code{QACollate} extends \code{.QA2}. It contains a \code{SimpleList}
  instance with zero or more \code{QASummary} elements.

  \code{QA} extends \code{.QA2}, and contains a \code{SimpleList} of
  zero or more \code{QASummary} elements. This class represents the
  results of the \code{qa2} analysis.

}

\section{Methods}{
  
  Methods defined on this class include:

  \describe{

    \item{qa2}{\code{signature(object="QACollate", state, ...,
	verbose=FALSE)} creates a QA report from the elements of
      \code{QACollate}. Methods on \code{qa2} for objects extending
      class \code{QASummary} summarize QA statistics for that class,
      e.g., \code{qa2,QAFrequentSequences-method} implements the
      calculations required to summarize frequently used sequences,
      using data in \code{state}.}

    \item{report}{\code{signature(x="QA", ...)} creates an HTML
      report. Methods on \code{report} for objects extending class
      \code{QASummary} are responsible for creating the html snippet for
      that QA component.}

    \item{flag}{\code{signature(object=".QA2", ..., verbose=FALSE)}
      implements criteria to flag individual lanes as failing quality
      assessment. NOTE: flag is not fully implemented.}

    \item{rbind}{\code{signature(...="QASummary")}: rbind multiple
      summary elements of the same class, as when these have been
      created by separately calculating statistics on a number of fastq
      files.}

    \item{show}{\code{signature(object = "SolexaExportQA")}: Display an
      overview of the object contents.}

  }
}

\seealso{
  Specific classes derived from \code{.QA2}
}

\author{Martin Morgan <mtmmorgan@fhcrc.org>}

\examples{
getClass(".QA2", where=getNamespace("ShortRead"))
}

\keyword{classes}
