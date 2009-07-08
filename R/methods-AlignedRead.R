###

.AlignedRead_validity <- function(object) {
    msg <- NULL
    len <- length(sread(object))
    slts <- c("chromosome", "position", "strand", "alignQuality")
    olen <- c(length(chromosome(object)),
              length(position(object)),
              length(strand(object)),
              length(alignQuality(object)))
    if (!all(olen==len)) {
        bad <- olen!=len
        msg <- c(msg,
                 sprintf("length mismatch: expected %d, found:\n  %s",
                         len, paste(slts[bad], olen[bad], sep="=",
                                    collapse=", ")))
    }
    if (is.null(msg)) TRUE else msg
}

setMethod(.srValidity, "AlignedRead", .AlignedRead_validity)

.AlignedRead_QualityConstructor <-
    function(sread)
{
    if (length(sread) > 0)
        unlist(lapply(width(sread), polyn, nucleotides="!"))
    else
        character(0)
}

AlignedRead <- function(sread = DNAStringSet(character(0)),
                        id = BStringSet(character(length(sread))),
                        quality = FastqQuality(
                          .AlignedRead_QualityConstructor(sread)),
                        chromosome = factor(rep(NA, length(sread))),
                        position = rep(NA_integer_, length(sread)),
                        strand = factor(rep(NA_integer_, length(sread)),
                          levels=.STRAND_LEVELS),
                        alignQuality = NumericQuality(
                          rep(NA_real_, length(sread))),
                        alignData = AlignedDataFrame(
                          nrow=length(sread)))
{
    new("AlignedRead", sread=sread, id=id, quality=quality,
        chromosome=as.factor(chromosome), position=position,
        strand=strand, alignQuality=alignQuality, alignData=alignData)
}

.make_getter(c("chromosome", "position", "alignQuality",
               "alignData"))

setMethod("strand", "AlignedRead", function(x)
{
    slot(x, "strand")
})

## coerce

setAs("PairwiseAlignedXStringSet", "AlignedRead",
      function(from, to) {
        pat <- pattern(from)
        quality <- character()
        if (is(pat, "QualityAlignedXStringSet"))
          quality <- quality(pat)
        new("AlignedRead", sread = unaligned(pat), id = names(pat),
            quality = FastqQuality(quality),
            position = start(Views(pat)),
            alignQuality = IntegerQuality(score(from)))
      })

## subset

setMethod("[", c("AlignedRead", "missing", "missing"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("AlignedRead", "missing", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("AlignedRead", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

.AlignedRead_subset <- function(x, i, j, ..., drop=TRUE) {
    if (nargs() != 2) .subset_err()
    initialize(x, sread=sread(x)[i], id=id(x)[i],
               quality=quality(x)[i],
               chromosome=factor(chromosome(x)[i]),
               position=position(x)[i], strand=strand(x)[i],
               alignQuality=alignQuality(x)[i],
               alignData=alignData(x)[i,]) }

setMethod("[", c("AlignedRead", "ANY", "missing"), .AlignedRead_subset)

setMethod(append, c("AlignedRead", "AlignedRead", "missing"),
    function(x, values, after=length(x))
{
    initialize(x,
               chromosome=.append.factor(chromosome(x),
                 chromosome(values)),
               position=append(position(x), position(values)),
               strand=.append.factor(strand(x), strand(values)),
               alignQuality=append(alignQuality(x),
                 alignQuality(values)),
               alignData=append(alignData(x), alignData(values)),
               quality=append(quality(x), quality(values)),
               sread=append(sread(x), sread(values)),
               id=append(id(x), id(values)))
})

## srorder, etc; srsort picked up by generic

setMethod(srorder, "AlignedRead", function(x, ...) {
    order(chromosome(x), strand(x), position(x), srorder(sread(x)))
})

setMethod(srrank, "AlignedRead", function(x, ...) {
    o <- srorder(x)
    .Call(.aligned_read_rank, x, o, environment())
})

setMethod(srduplicated, "AlignedRead", function(x, ...) {
    duplicated(srrank(x, ...))
})

## coverage

setMethod(coverage, "AlignedRead",
    function(x, start=NA, end=NA, shift=0L, width=NULL, weight=1L, ...,
             coords=c("leftmost", "fiveprime"),
             extend=0L)
{
##---- remove this block when the start/end interface gets finally dropped ---
    isSingleNA <- function(x) {is.atomic(x) && length(x) == 1 && is.na(x)}
    if (!isSingleNA(start) || !isSingleNA(end)) {
        if (!identical(shift, 0L) || !is.null(width)) {
            msg <- "use only 'shift'/'width' (or 'start'/'end')"
            .throw(SRError("UserArgumentMismatch", msg))
        }
        msg <- "use 'shift'/'width'  instead of 'start'/'end'"
        .throw(SRWarn("UserArgumentMismatch", msg))
        if (!is.numeric(start) || any(is.na(start))
         || !is.numeric(end) || any(is.na(end))) {
            msg <- "'start', 'end' must be named integer values"
            .throw(SRError("UserArgumentMismatch", msg))
        }
        if (!all(names(start) == names(end)))
            .throw(SRError("UserArgumentMismatch", 
                           "not all names(start) == names(end)"))
        shift <- 1L - start
        width <- end + shift
        names(shift) <- names(width) <- names(start)
        if (any(width < 0L)) {
            .throw(SRError("UserArgumentMismatch", 
                           "'end' must be >= 'start' - 1"))
        }
    }
##----------------------------------------------------------------------------

    ## Argument checking:
    chrlvls <- levels(chromosome(x))
    if (!identical(shift, 0L)) {
        if (!is.numeric(shift)) {
            .throw(SRError("UserArgumentMismatch",
                           "if '%s' is not 0L, then it must be a vector of integer values\n  see %s",
                           "shift", '?"AlignedRead-class"'))
        }
        if (!all(chrlvls %in% names(shift))) { 
            .throw(SRError("UserArgumentMismatch",
                           "'names(%s)' (or 'names(%s)') mismatch with 'levels(chromosome(x))'\n  see %s",
                           "shift", "start", '?"AlignedRead-class"'))
        }
        if (any(duplicated(names(shift)))) {
            .throw(SRError("UserArgumentMismatch",
                           "'names(%s)' (or 'names(%s)') have duplicates\n  see %s",
                           "shift", "start", '?"AlignedRead-class"'))
        }
    }
    if (!is.null(width)) {
        if (!is.numeric(width)) {
            .throw(SRError("UserArgumentMismatch",
                           "if '%s' is not NULL, then it must be a vector of integer values\n  see %s",
                           "width", '?"AlignedRead-class"'))
        }
        if (!all(chrlvls %in% names(width))) {
            .throw(SRError("UserArgumentMismatch",
                           "'names(%s)' (or 'names(%s)') mismatch with 'levels(chromosome(x))'\n  see %s",
                           "width", "end", '?"AlignedRead-class"'))
        }
        if (any(duplicated(names(width)))) {
            .throw(SRError("UserArgumentMismatch",
                           "'names(%s)' (or 'names(%s)') have duplicates\n  see %s",
                           "width", "end", '?"AlignedRead-class"'))
        }
    }
    if (!identical(weight, 1L)) {
        .throw(SRError("UserArgumentMismatch",
                       "weighting the reads is not supported yet, sorry\n  see %s",
                       '?"AlignedRead-class"'))
    }
    tryCatch(coords <- match.arg(coords),
        error=function(err) {
            vals <- formals(sys.function(sys.parent(4)))[["coords"]]
            .throw(SRError("UserArgumentMismatch",
                           "'%s' must be one of '%s'\n  see %s", "coords",
                           paste(eval(vals), collapse="' '"),
                           '?"AlignedRead-class"'))
        })
    if (!is.integer(extend) ||
        !(length(extend) == 1 || length(extend) == length(x))) {
        .throw(SRError("UserArgumentMismatch",
                       "'%s' must be '%s'", "extend",
                       "integer(n)', n == 1 or length(x)"))
    }
    ## end of argument checking.

    if (coords == "leftmost") {
        rstart <- position(x) -
            ifelse(strand(x) == "+", 0L, extend)
        rend <- position(x) + width(x) - 1L + 
            ifelse(strand(x) == "+", extend, 0L)
    } else {
        rstart <- position(x) -
            ifelse(strand(x) == "+", 0L, width(x) + extend - 1L)
        rend <- position(x) +
            ifelse(strand(x) == "+", width(x) + extend - 1L, 0L)
    }
    cvg <- lapply(chrlvls,
                  function(chr, ...) {
                      idx <- chromosome(x) == chr
                      chr_rstart <- rstart[idx]
                      chr_rend <- rend[idx]
                      if (identical(shift, 0L))
                          chr_shift <- 0L
                      else
                          chr_shift <- shift[chr]
                      if (is.null(width))
                          chr_width <- NULL
                      else
                          chr_width <- width[chr]
                      coverage(IRanges(chr_rstart, chr_rend),
                               shift=chr_shift, width=chr_width, ...)
                  },
                  ...
           )
    names(cvg) <- chrlvls
    GenomeData(cvg, method="coverage,AlignedRead-method",
               coords=coords, extend=extend)
})

## show

setMethod(show, "AlignedRead", function(object) {
    callNextMethod()
    cat("chromosome:", selectSome(chromosome(object)), "\n")
    cat("position:", selectSome(position(object)), "\n")
    cat("strand:", selectSome(strand(object)), "\n")
    cat("alignQuality:", class(alignQuality(object)), "\n")
    cat("alignData varLabels:",
        selectSome(varLabels(alignData(object))), "\n")
})

setMethod(detail, "AlignedRead", function(object, ...) {
    callNextMethod()
    cat("\nchromosome:", selectSome(chromosome(object)), "\n")
    cat("position:", selectSome(position(object)), "\n")
    cat("strand:", selectSome(strand(object)), "\n")
    cat("alignQuality:\n")
    detail(alignQuality(object))
    cat("\nalignData:\n")
    show(alignData(object))
})

## summary

## perhaps summary statistics like ShortReadQ except broken down by chromosome,
## strand, and their combination
