.AlignedRead_validity <- function(object) {
    msg <- NULL
    len <- length(sread(object))
    slts <- c("chromosome", "position", "strand", "alignQuality")
    olen <- sapply(slts, function(elt) {
        length(do.call(elt, list(object)))
    })
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

setMethod(strand, "AlignedRead", function(x, ...) 
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
    tryCatch(coords <- match.arg(coords),
        error=function(err) {
            vals <- formals(sys.function(sys.parent(4)))[["coords"]]
            .throw(SRError("UserArgumentMismatch",
                           "'%s' must be one of '%s'\n  see %s", "coords",
                           paste(eval(vals), collapse="' '"),
                           '?"AlignedRead-class"'))
        })
    chrlvls <- levels(chromosome(x))
    if (length(start) != 1 && 
        !all(c(names(start) %in% chrlvls,
               chrlvls %in% names(start)))) {
        .throw(SRError("UserArgumentMismatch",
                       "'names(%s)' mismatch with 'levels(chromosome(x))'\n  see %s",
                       "start", '?"AlignedRead-class"'))
    }
    if (length(end) != 1 &&
        !all(c(names(end) %in% chrlvls,
               chrlvls %in% names(end)))) {
        .throw(SRError("UserArgumentMismatch",
                       "'names(%s)' mismatch with 'levels(chromosome(x))'\n  see %s",
                       "end",  '?"AlignedRead-class"'))
    }
    if (!is.integer(extend) ||
        !(length(extend) == 1 || length(extend) == length(x))) {
        .throw(SRError("UserArgumentMismatch",
                       "'%s' must be '%s'", "extend",
                       "integer(n)', n == 1 or length(x)"))
    }
    if (coords == "leftmost") {
        rstart <- ifelse(strand(x)=="+", position(x),
                         position(x) - extend)
        rend <- ifelse(strand(x) == "+",
                       position(x) + width(x) + extend - 1L,
                       position(x) + width(x) - 1L)
    } else {
        rstart <- ifelse(strand(x) == "+", position(x),
                         position(x) - width(x) - extend + 1L)
        rend <- ifelse(strand(x) == "+",
                       position(x) + width(x) + extend - 1L,
                       position(x))
    }
    cvg <- lapply(chrlvls, function(chr, aln, rstart, rend, start, end, ...) {
        idx <- chromosome(aln) == chr
        rstart <- rstart[idx]
        rend <- rend[idx]
        if (length(start) == 1 && is.na(start)) start <- min(rstart)
        else if (length(start) != 1) start <- start[chr]
        if (length(end) == 1 && is.na(end)) end <- max(rend)
        else if (length(end) != 1) end <- end[chr]
        coverage(IRanges(rstart, rend),
                 shift=1L-start, width=end+1L-start, ...)
    }, x, rstart, rend, start, end, ...)
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
