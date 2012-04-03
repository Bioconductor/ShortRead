\name{qa2}

\alias{QAFastqSource}

\alias{QACollate}
\alias{QA}

\alias{QAFlagged}
\alias{QAFiltered}
\alias{QAAdapterContamination}
\alias{QAData}
\alias{QAFrequentSequence}
\alias{QANucleotideByCycle}
\alias{QANucleotideUse}
\alias{QAQualityByCycle}
\alias{QAQualityUse}
\alias{QAReadQuality}
\alias{QASequenceUse}

\alias{QACollate,QAFastqSource-method}
\alias{QACollate,missing-method}

\alias{qa2}
\alias{qa2,FastqSampler-method}
\alias{qa2,QAAdapterContamination-method}
\alias{qa2,QACollate-method}
\alias{qa2,QAFastqSource-method}
\alias{qa2,QAFrequentSequence-method}
\alias{qa2,QANucleotideByCycle-method}
\alias{qa2,QANucleotideUse-method}
\alias{qa2,QAQualityByCycle-method}
\alias{qa2,QAQualityUse-method}
\alias{qa2,QAReadQuality-method}
\alias{qa2,QASequenceUse-method}

\alias{flag}
\alias{flag,.QA2-method}
\alias{flag,QAFrequentSequence-method}
\alias{flag,QAReadQuality-method}
\alias{flag,QASource-method}

\alias{report,QA-method}
\alias{report,QAAdapterContamination-method}
\alias{report,QAFiltered-method}
\alias{report,QAFlagged-method}
\alias{report,QAFrequentSequence-method}
\alias{report,QANucleotideByCycle-method}
\alias{report,QANucleotideUse-method}
\alias{report,QAQualityByCycle-method}
\alias{report,QAQualityUse-method}
\alias{report,QAReadQuality-method}
\alias{report,QASequenceUse-method}
\alias{report,QASource-method}

\alias{rbind,QASummary-method}

\alias{show,QAAdapterContamination-method}
\alias{show,QACollate-method}
\alias{show,QAFastqSource-method}
\alias{show,QAFrequentSequence-method}
\alias{show,QAReadQuality-method}
\alias{show,QASummary-method}

\title{(Updated) quality assessment reports on short reads}

\description{

  This page summarizes an updated approach to quality assessment reports
  in \code{ShortRead}.

}
\usage{
## Input source for short reads
QAFastqSource(con = character(), n = 1e+06, readerBlockSize = 1e+08,
    flagNSequencesRange = NA_integer_, ...,
    html = system.file("template", "QASources.html", package="ShortRead"))
QAData(seq = ShortReadQ(), filter = logical(length(seq)), ...)

## Possible QA elements
QAFrequentSequence(useFilter = TRUE, addFilter = TRUE,
    n = NA_integer_, a = NA_integer_, flagK=.8, reportSequences = FALSE,
    ...)
QANucleotideByCycle(useFilter = TRUE, addFilter = TRUE, ...)
QANucleotideUse(useFilter = TRUE, addFilter = TRUE, ...)
QAQualityByCycle(useFilter = TRUE, addFilter = TRUE, ...)
QAQualityUse(useFilter = TRUE, addFilter = TRUE, ...)
QAReadQuality(useFilter = TRUE, addFilter = TRUE,
    flagK = 0.2, flagA = 30L, ...)
QASequenceUse(useFilter = TRUE, addFilter = TRUE, ...)
QAAdapterContamination(useFilter = TRUE, addFilter = TRUE,
    Lpattern = NA_character_, Rpattern = NA_character_,
    max.Lmismatch = 0.1, max.Rmismatch = 0.2, min.trim = 9L, ...)

## Order QA report elements
QACollate(src, ...)

## perform analysis
qa2(object, state, ..., verbose=FALSE)

## Outputs from qa2
QA(src, filtered, flagged, ...)
QAFiltered(useFilter = TRUE, addFilter = TRUE, ...)
QAFlagged(useFilter = TRUE, addFilter = TRUE, ...)

## Summarize results as html report
\S4method{report}{QA}(x, ..., dest = tempfile(), type = "html")

## additional methods; 'flag' is not fully implemented
flag(object, ..., verbose=FALSE)

\S4method{rbind}{QASummary}(..., deparse.level = 1)
}

\arguments{

  \item{con}{\code{character(1)} file location of fastq input, as used
    by \code{FastqSampler}.}

  \item{n}{\code{integer(1)} number of records to input, as used by
    \code{FastqStreamer} (\code{QAFastqSource}). \code{integer(1)}
    number of sequences to tag as \sQuote{frequent}
    (\code{QAFrequentSequence}). }

  \item{readerBlockSize}{integer(1) number of bytes to input, as used by
    \code{FastqStreamer}.}

  \item{flagNSequencesRange}{\code{integer(2)} minimum and maximum reads
    above which source files will be flagged as outliers.}

  \item{html}{\code{character(1)} location of the HTML template for
    summarizing this report element.}

  \item{seq}{\code{\link{ShortReadQ}} representation of fastq data.}

  \item{filter}{\code{logical()} vector with length equal to
    \code{seq}, indicating whether elements of \code{seq} are filtered
    (\code{TRUE}) or not.}

  \item{useFilter, addFilter}{\code{logical(1)} indicating whether the
    QA element should be calculating using the filtered
    (\code{useFilter=TRUE}) or all reads, and whether reads failing the
    QA element should be added to the filter used by subsequent steps
    (\code{addFilter = TRUE}) or not.}

  \item{a}{\code{integer(1)} count of number of sequences above which a
    read will be considered \sQuote{frequent}
    (\code{QAFrequentSequence}).}

  \item{flagK, flagA}{\code{flagK} \code{numeric(1)} between 0 and 1
    indicating the fraction of frequent sequences greater than or
    equal to \code{n} or \code{a} above which a fastq file will be
    flagged (\code{QAFrequentSequence}). \code{flagK} \code{numeric{1}}
    between 0 and 1 and \code{flagA} \code{integer(1)} indicating that a
    run should be flagged when the fraction of reads with quality
    greater than or equal to \code{flagA} falls below threshold
    \code{flagK}.}

  \item{reportSequences}{\code{logical(1)} indicating whether frequent
    sequences are to be reported.}

  \item{Lpattern, Rpattern, max.Lmismatch, max.Rmismatch,
    min.trim}{Parameters influencing adapter identification, see
    \code{\link{matchPattern}}.}

  \item{src}{The source, e.g., \code{QAFastqSource}, on which the
    quality assessment report will be based.}

  \item{object}{An instance of class derived from \code{QA} on which
    quality metrics will be derived; for end users, this is usually the
    result of \code{QACollate}.}.

  \item{state}{The data on which quality assessment will be performed;
    this is not usually necessary for end-users.}

  \item{verbose}{\code{logical(1)} indicating whether progress reports
    should be reported.}

  \item{filtered, flagged}{Primarily for internal use, instances of
    \code{QAFiltered} and \code{QAFlagged}.}

  \item{x}{An instance of \code{QA} on which a report is to be
    generated.}

  \item{dest}{\code{character(1)} providing the directory in which the
    report is to be generated.}

  \item{type}{\code{character(1)} indicating the type of report to be
    generated; only \dQuote{html} is supported.}

  \item{deparse.level}{see \code{\link{rbind}}.}

  \item{...}{Additional arguments, e.g., \code{html} to specify the
    location of the html source to use as a template for the report.}

}

\details{

  Use \code{QACollate} to specify an order in which components of a QA
  report are to be assembled. The first argument is the data source
  (e.g., \code{QAFastqSource}).

  Functions related to data input include:
  \describe{

    \item{\code{QAFastqSource}}{defines the location of fastq files to
      be included in the report. \code{con} is used to construct a
      \code{\link{FastqSampler}} instance, and records are processed
      using \code{qa2,QAFastqSource-method}.}

    \item{\code{QAData}}{is a class for representing the data during the
      QA report generation pass; it is primarily for internal use.}

  }

  Possible elements in a QA report are:
  \describe{

    \item{\code{QAFrequentSequence}}{identifies the most-commonly
      occuring sequences. One of \code{n} or \code{a} can be non-NA, and
      determine the number of frequent sequences reported.  \code{n}
      specifies the number of most-frequent sequences to filter, e.g.,
      \code{n=10} would filter the top 10 most commonly occurring
      sequences; \code{a} provides a threshold frequency (count) above
      which reads are filtered. The sample is flagged when a fraction
      \code{flagK} of the reads are filtered.

      \code{reportSequences} determines whether the most commonly
      occuring sequences, as determined by \code{n} or \code{a}, are
      printed in the html report.
    }

    \item{\code{QANucleotideByCycle}}{reports nucleotide frequency as a
      function of cycle.}

    \item{\code{QAQualityByCycle}}{reports average quality score as a
      function of cycle.}

    \item{\code{QAQualityUse}}{summarizes overall nucleotide qualities.}

    \item{\code{QAReadQuality}}{summarizes the distribution of read
      qualities.}

    \item{\code{QASequenceUse}}{summarizes the cumulative distribution
      of reads occurring 1, 2, \dots times.}

    \item{\code{QAAdapterContamination}}{reports the occurrence of
      \sQuote{adapter} sequences on the left and / or right end of each
      read.}

  }
}

\value{

  An object derived from class \code{\linkS4class{.QA}}. Values
  contained in this object are meant for use by \code{\link{report}}

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{\code{\linkS4class{QA}}.}

\examples{
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")
fls <- c(fl, fl)

coll <- QACollate(QAFastqSource(fls), QAReadQuality(),
    QAAdapterContamination(), QANucleotideUse(),
    QAQualityUse(), QASequenceUse(),
    QAFrequentSequence(n=10), QANucleotideByCycle(),
    QAQualityByCycle())
x <- qa2(coll, verbose=TRUE)
res <- report(x)
if (interactive()) 
    browseURL(res)

\dontrun{## parallel evaluation
options(srapply_fapply="parallel")
coll <- QACollate(QAFastqSource(fls), qaFastqTemplate)
x2 <- qa2(coll, verbose=TRUE)
browseURL(res2 <- report(x2))
}

}

\keyword{manip}
