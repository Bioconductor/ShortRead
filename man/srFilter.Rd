\name{srFilter}
\alias{srFilter}
\alias{srFilter,missing-method}
\alias{srFilter,function-method}
\alias{idFilter}
\alias{chromosomeFilter}
\alias{positionFilter}
\alias{strandFilter}
\alias{occurrenceFilter}
\alias{nFilter}
\alias{polynFilter}
\alias{dustyFilter}
\alias{srdistanceFilter}
\alias{alignQualityFilter}
\alias{alignDataFilter}
\alias{compose}

\title{Functions for user-created and built-in ShortRead filters}

\description{

  These functions create user-defined (\code{srFitler}) or built-in
  instances of \code{\linkS4class{SRFilter}}  objects. Filters can be
  applied to objects from \code{ShortRead}, returning a logical vector
  to be used to subset the objects to include only those components
  satisfying the filter.

}
\usage{
srFilter(fun, name = NA_character_, ...)
\S4method{srFilter}{missing}(fun, name=NA_character_, ...)
\S4method{srFilter}{function}(fun, name=NA_character_, ...)

compose(filt, ..., .name)

idFilter(regex=character(0), fixed=FALSE, exclude=FALSE,
         .name="idFilter")
occurrenceFilter(min=1L, max=1L,
                 withSread=c(NA, TRUE, FALSE),
                 duplicates=c("head", "tail", "sample", "none"),
                 .name=.occurrenceName(min, max, withSread,
                                       duplicates))
nFilter(threshold=0L, .name="CleanNFilter")
polynFilter(threshold=0L, nuc=c("A", "C", "T", "G", "other"),
           .name="PolyNFilter")
dustyFilter(threshold=Inf, batchSize=NA, .name="DustyFilter")
srdistanceFilter(subject=character(0), threshold=0L,
                 .name="SRDistanceFilter")

##
## legacy filters for ungapped alignments
##

chromosomeFilter(regex=character(0), fixed=FALSE, exclude=FALSE,
                 .name="ChromosomeFilter")
positionFilter(min=-Inf, max=Inf, .name="PositionFilter")
strandFilter(strandLevels=character(0), .name="StrandFilter")
alignQualityFilter(threshold=0L, .name="AlignQualityFilter")
alignDataFilter(expr=expression(), .name="AlignDataFilter")
}

\arguments{

  \item{fun}{An object of class \code{function} to be used as a
    filter. \code{fun} must accept a single named argument \code{x}, and
    is expected to return a logical vector such that \code{x[fun(x)]}
    selects only those elements of \code{x} satisfying the conditions of
    \code{fun}
  }

  \item{name}{A \code{character(1)} object to be used as the name of the
    filter.  The \code{name} is useful for debugging and reference.}

  \item{filt}{A \code{\linkS4class{SRFilter}} object, to be used with
    additional arguments to create a composite filter.}

  \item{.name}{An optional \code{character(1)} object used to over-ride
    the name applied to default filters.}

  \item{regex}{Either \code{character(0)} or a \code{character(1)}
    regular expression used as \code{grep(regex, chromosome(x))} to
    filter based on chromosome. The default (\code{character(0)})
    performs no filtering}

  \item{fixed}{\code{logical(1)} passed to \code{\link{grep}},
    influencing how pattern matching occurs.}

  \item{exclude}{\code{logical(1)} which, when \code{TRUE}, uses
    \code{regex} to exclude, rather than include, reads.}

  \item{min}{\code{numeric(1)}}

  \item{max}{\code{numeric(1)}. For \code{positionFilter}, \code{min}
    and \code{max} define the closed interval in which position must be
    found \code{min <= position <= max}. For \code{occurrenceFilter},
    \code{min} and \code{max} define the minimum and maximum number of
    times a read occurs after the filter.}

  \item{strandLevels}{Either \code{character(0)} or \code{character(1)}
    containing strand levels to be selected. \code{ShortRead} objects
    have standard strand levels \code{NA, "+", "-", "*"}, with \code{NA}
    meaning strand information not available and \code{"*"} meaning
    strand information not relevant.}

  \item{withSread}{A \code{logical(1)} indicating whether uniqueness
    includes the read sequence (\code{withSread=TRUE}), is based only on
    chromosome, position, and strand (\code{withSread=FALSE}), or only
    the read sequence (\code{withSread=NA}), as described for
    \code{occurrenceFilter} below..}

  \item{duplicates}{Either \code{character{1}}, a function \code{name},
	or a function taking a single argument. Influence how duplicates are
	handled, as described for \code{occurrenceFilter} below.}

  \item{threshold}{A \code{numeric(1)} value representing a minimum
    (\code{srdistanceFilter}, \code{alignQualityFilter}) or maximum
    (\code{nFilter}, \code{polynFilter}, \code{dustyFilter}) criterion
    for the filter. The minima and maxima are closed-interval (i.e.,
    \code{x >= threshold}, \code{x <= threshold} for some property
    \code{x} of the object being filtered).}

  \item{nuc}{A \code{character} vector containing IUPAC symbols for
    nucleotides or the value \code{"other"} corresponding to all
    non-nucleotide symbols, e.g., \code{N}.}

  \item{batchSize}{\code{NA} or an \code{integer(1)} vector indicating
    the number of DNA sequences to be processed simultaneously by
    \code{dustyFilter}. By default, all reads are processed
    simultaneously. Smaller values use less memory but are
    computationally less efficient.}

  \item{subject}{A \code{character()} of any length, to be used as the
    corresponding argument to \code{\link{srdistance}}.}

  \item{expr}{A \code{expression} to be evaluated with
    \code{pData(alignData(x))}.}

  \item{\dots}{Additional arguments for subsequent methods; these
    arguments are not currently used.}
}

\details{

  \code{srFilter} allows users to construct their own filters. The
  \code{fun} argument to \code{srFilter} must be a function accepting a
  single argument \code{x} and returning a logical vector that can be
  used to select elements of \code{x} satisfying the filter with
  \code{x[fun(x)]}

  The \code{signature(fun="missing")} method creates a default filter
  that returns a vector of \code{TRUE} values with length equal to
  \code{length(x)}.

  \code{compose} constructs a new filter from one or more existing
  filter. The result is a filter that returns a logical vector with
  indices corresponding to components of \code{x} that pass all
  filters. If not provided, the name of the filter consists of the names
  of all component filters, each separated by \code{" o "}.

  The remaining functions documented on this page are built-in filters
  that accept an argument \code{x} and return a logical vector of
  \code{length(x)} indicating which components of \code{x} satisfy the
  filter.

  \code{idFilter} selects elements satisfying
  \code{grep(regex, id(x), fixed=fixed)}.

  \code{chromosomeFilter} selects elements satisfying
  \code{grep(regex, chromosome(x), fixed=fixed)}.

  \code{positionFilter} selects elements satisfying
  \code{min <= position(x) <= max}.

  \code{strandFilter} selects elements satisfying
  \code{match(strand(x), strand, nomatch=0) > 0}.

  \code{occurrenceFilter} selects elements that occur \code{>=min} and
  \code{<=max} times. \code{withSread} determines how reads will be
  treated: \code{TRUE} to include the sread, chromosome, strand, and
  position when determining occurrence, \code{FALSE} to include
  chromosome, strand, and position, and \code{NA} to include only
  sread. The default is \code{withSread=NA}. \code{duplicates}
  determines how reads with more than \code{max} reads are
  treated. \code{head} selects the first \code{max} reads of each set of
  duplicates, \code{tail} the last \code{max} reads, and \code{sample} a
  random sample of \code{max} reads. \code{none} removes all reads
  represented more than \code{max} times. The user can also provide a
  function (as used by \code{\link{tapply}}) of a single argument to
  select amongst reads.
  

  \code{nFilter} selects elements with fewer than \code{threshold}
  \code{'N'} symbols in each element of \code{sread(x)}.

  \code{polynFilter} selects elements with fewer than \code{threshold}
  copies of any nucleotide indicated by \code{nuc}.

  \code{dustyFilter} selects elements with high sequence complexity, as
  characterized by their \code{\link{dustyScore}}. This emulates the
  \code{dust} command from \code{WindowMaker}
  software. Calculations can be memory intensive; use
  \code{batchSize} to process the argument to \code{dustyFilter} in
  batches of the specified size.
  
  \code{srdistanceFilter} selects elements at an edit distance greater
  than \code{threshold} from all sequences in \code{subject}.

  \code{alignQualityFilter} selects elements with \code{alignQuality(x)}
  greater than \code{threshold}.

  \code{alignDataFilter} selects elements with
  \code{pData(alignData(x))} satisfying \code{expr}. \code{expr} should
  be formulated as though it were to be evaluated as
  \code{eval(expr, pData(alignData(x)))}.

}

\value{

  \code{srFilter} returns an object of \code{\linkS4class{SRFilter}}.

  Built-in filters return a logical vector of \code{length(x)}, with
  \code{TRUE} indicating components that pass the filter.

}

\author{Martin Morgan <mtmorgan@fhcrc.org>}

\seealso{\code{\linkS4class{SRFilter}}.}

\examples{
sp <- SolexaPath(system.file("extdata", package="ShortRead"))
aln <- readAligned(sp, "s_2_export.txt") # Solexa export file, as example

# a 'chromosome 5' filter
filt <- chromosomeFilter("chr5.fa")
aln[filt(aln)]
# filter during input
readAligned(sp, "s_2_export.txt", filter=filt)

# x- and y- coordinates stored in alignData, when source is SolexaExport
xy <- alignDataFilter(expression(abs(x-500) > 200 & abs(y-500) > 200))
aln[xy(aln)]

# both filters as a single filter
chr5xy <- compose(filt, xy)
aln[chr5xy(aln)]

# both filters as a collection
filters <- c(filt, xy)
subsetByFilter(aln, filters)
summary(filters, aln)

# read, chromosome, strand, position tuples occurring exactly once
aln[occurrenceFilter(withSread=TRUE, duplicates="none")(aln)]
# reads occurring exactly once
aln[occurrenceFilter(withSread=NA, duplicates="none")(aln)]
# chromosome, strand, position tuples occurring exactly once
aln[occurrenceFilter(withSread=FALSE, duplicates="none")(aln)]

# custom filter: minimum calibrated base call quality >20
goodq <- srFilter(function(x) {
    apply(as(quality(x), "matrix"), 1, min, na.rm=TRUE) > 20
}, name="GoodQualityBases")
goodq
aln[goodq(aln)]
}

\keyword{manip}
