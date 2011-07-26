### =========================================================================
### GappedReads objects
### -------------------------------------------------------------------------
###

setClass("GappedReads",
    contains="GappedAlignments",
    representation(
        qseq="DNAStringSet"
        ## TODO: Maybe add the read quality? mismatch information?
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setGeneric("qseq", function(x) standardGeneric("qseq"))

setMethod("qseq", "GappedReads", function(x) x@qseq)

### Overriding "qwidth" method for GappedAlignments objects with a faster
### method.
setMethod("qwidth", "GappedReads", function(x) width(qseq(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GappedReads.qseq <- function(x)
{
    x_qseq <- qseq(x)
    if (class(x_qseq) != "DNAStringSet" || !is.null(names(x_qseq)))
        return("'qseq(x)' must be an unnamed DNAStringSet instance")
    if (length(x_qseq) != length(cigar(x)))
        return("'qseq(x)' and 'cigar(x)' must have the same length")
    if (!identical(width(x_qseq), cigarToQWidth(cigar(x))))
        return(paste("'width(qseq(x))' and 'cigarToQWidth(cigar(x))'",
                     "must be identical"))
    NULL
}

.valid.GappedReads <- function(x)
{
    .valid.GappedReads.qseq(x)
}

setValidity2("GappedReads", .valid.GappedReads,
             where=asNamespace("ShortRead"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

GappedReads <- function(rname=Rle(factor()), pos=integer(0),
                        cigar=character(0), strand=NULL,
                        qseq=DNAStringSet(),
                        names=NULL, seqlengths=NULL)
{
    galn <- GappedAlignments(rname=rname, pos=pos,
                             cigar=cigar, strand=strand,
                             names=names, seqlengths=seqlengths)
    new("GappedReads", galn, qseq=qseq)
}

readGappedReads <- function(file, format="BAM", use.names=FALSE, ...)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if (!isSingleString(format))
        stop("'format' must be a single string")
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    dotargs <- list(...)
    if (length(dotargs) != 0L && is.null(names(dotargs)))
        stop("extra arguments must be named")
    if (format == "BAM") {
        if ("index" %in% names(dotargs)) {
            index <- dotargs$index
            dotargs$index <- NULL
        } else {
            index <- file
        }
        if ("which" %in% names(dotargs)) {
            which <- dotargs$which
            dotargs$which <- NULL
        } else {
            which <- RangesList()
        }
        args <- c(list(file=file, index=index, use.names=use.names),
                  dotargs,
                  list(which=which))
        ans <- do.call(readBamGappedReads, args)
        return(ans)
    }
    stop("only BAM format is supported for now")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "GappedReads",
    function(x, i, j, ... , drop=TRUE)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (is(i, "Rle"))
            i <- as.vector(i)
        if (!is.atomic(i))
            stop("invalid subscript type")
        lx <- length(x)
        if (length(i) == 0L) {
            i <- integer(0)
        } else if (is.numeric(i)) {
            if (min(i) < 0L)
                i <- seq_len(lx)[i]
            else if (!is.integer(i))
                i <- as.integer(i)
        } else if (is.logical(i)) {
            if (length(i) > lx)
                stop("subscript out of bounds")
            i <- seq_len(lx)[i]
        } else {
            stop("invalid subscript type")
        }
        x <- callNextMethod()
        x@qseq <- x@qseq[i]
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining.
###

setMethod("c", "GappedReads",
    function (x, ..., recursive = FALSE)
    {
        stop("coming soon")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "qnarrow", "narrow", and "pintersect" methods.
###

setMethod("qnarrow", "GappedReads",
    function(x, start=NA, end=NA, width=NA)
    {
        stop("coming soon")
        ans_cigar <- cigarQNarrow(cigar(x),
                                  start=start, end=end, width=width)
        ans_start <- start(x) + attr(ans_cigar, "rshift")
        updateCigarAndStart(x, cigar=ans_cigar, start=ans_start)
    }
)

setMethod("narrow", "GappedReads",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        stop("coming soon")
        ans_cigar <- cigarNarrow(cigar(x),
                                 start=start, end=end, width=width)
        ans_start <- start(x) + attr(ans_cigar, "rshift")
        updateCigarAndStart(x, cigar=ans_cigar, start=ans_start)
    }
)

