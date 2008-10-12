\name{srFilter}
\alias{srFilter}
\alias{srFilter,missing-method}
\alias{srFilter,function-method}
\alias{chromosomeFilter}
\alias{strandFilter}
\alias{nFilter}
\alias{polynFilter}
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
\S4method{srFilter}{missing}(fun, name=NA_character_)
\S4method{srFilter}{function}(fun, name=NA_character_)

compose(filt, ..., .name)

chromosomeFilter(regex=character(0), .name="ChromosomeFilter")
strandFilter(strandLevels=character(0), .name="StrandFilter")
nFilter(threshold=0L, .name="CleanNFilter")
polynFilter(threshold=0L, nuc=c("A", "C", "T", "G", "other"),
           .name="PolyNFilter")
srdistanceFilter(subject=character(0), threshold=0L,
                 .name="SRDistanceFilter")
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
    additional arugments to create a composite filter.}

  \item{.name}{An optional \code{character(1)} object used to over-ride
    the name applied to default filters.}

  \item{regex}{Either \code{character(0)} or a \code{character(1)}
    regular expression used as \code{grep(regex, chromosome(x))} to
    filter based on chromosome. The default (\code{character(0)})
    performs no filtering}

  \item{strandLevels}{Either \code{character(0)} or \code{character(1)}
    containing strand levels to be selected. \code{ShortRead} objects
    have standard strand levels \code{NA, "+", "-", "*"}, with \code{NA}
    meaning strand information not available and \code{"*"} meaning
    strand information not relevant.}

  \item{threshold}{A \code{numeric(1)} value representing a minimum
    (\code{srdistanceFilter}, \code{alignQualityFilter}) or maximum
    (\code{nFilter}, \code{polynFilter}) criterion for the filter. The
    minima and maxima are closed-interval (i.e., \code{x >= threshold},
    \code{x <= threshold} for some property \code{x} of the object being
    filtered).}

  \item{nuc}{A \code{character} vector containing IUPAC symbols for
    nucleotides or the value \code{"other"} corresponding to all
    non-nucleotide symbols, e.g., \code{N}.}

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
  indicies corresponding to components of \code{x} that pass all
  filters. If not provided, the name of the filter consists of the names
  of all component filters, each separated by \code{" o "}.

  The remaining functions documented on this page are built-in filters
  that accept an argument \code{x} and return a logical vector of
  \code{length(x)} indicating which components of \code{x} satisfy the
  filter.

  \code{chromosomeFilter} selects elements satisfying
  \code{grep(regex, chromosome(x))}.

  \code{strandFilter} selects elemenst satisfying
  \code{match(strand(x), strand, nomatch=0) > 0}.

  \code{nFilter} selects elements with fewer than \code{threshold}
  \code{'N'} symbols in each element of \code{sread(x)}.

  \code{polynFilter} selects elements with fewer than \code{threshold}
  copies of any nucleotide indicated by \code{nuc}.

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
chr5 <- chromosomeFilter("chr5.fa")
aln[chr5(aln)]

# x- and y- coordinates stored in alignData, when source is SolexaExport
xy <- alignDataFilter(expression(abs(x-500) > 200 & abs(y-500) > 200))
aln[xy(aln)]

# both filters
chr5xy <- compose(chr5, xy)
aln[chr5xy(aln)]

# custom filter: minimum calibrated base call quality >20
goodq <- srFilter(function(x) {
    apply(as(quality(x), "matrix"), 1, min) > 20
}, name="GoodQualityBases")
goodq
aln[goodq(aln)]
}

\keyword{manip}
