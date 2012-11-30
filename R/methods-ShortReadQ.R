## validity / accessors / constructors

setMethod(.srValidity, "ShortReadQ", function(object) {
    msg <- NULL
    lenq <- length(quality(object))
    lens <- length(sread(object))
    if (lenq != lens) {
        txt <- sprintf("sread and quality length mismatch: %d %d",
                       lenq, lens)
        msg <- c(msg, txt)
    }
    if (!all(width(quality(object)) == width(sread(object)))) {
        txt <- sprintf("some sread and quality widths differ")
        msg <- c(msg, txt)
    }
    if (is.null(msg)) TRUE else msg
})

setMethod(ShortReadQ, c("DNAStringSet", "QualityScore", "BStringSet"),
          function(sread, quality, id, ...)
{
    new("ShortReadQ", sread=sread, quality=quality, id=id, ...)
})

setMethod(ShortReadQ, c("DNAStringSet", "QualityScore", "missing"),
          function(sread, quality, id, ...)
{
    new("ShortReadQ", sread=sread, quality=quality,
        id=BStringSet(character(length(sread))), ...)
})

setMethod(ShortReadQ, c("DNAStringSet", "BStringSet", "BStringSet"),
    function(sread, quality, id, ...,
             qualityType=c("Auto", "FastqQuality", "SFastqQuality"),
             filter=srFilter(), withIds=TRUE)
{
    if (!missing(filter))
        .check_type_and_length(filter, "SRFilter", NA)
    tryCatch({
        qualityType <- match.arg(qualityType)
    }, error=function(err) {
        .throw(SRError("UserArgumentMismatch", conditionMessage(err)))
    })
    tryCatch({
        qualityFunc <-
            switch(qualityType, Auto={
                alf <- alphabetFrequency(head(quality, 10000),
                                         collapse=TRUE)
                if (any(alf) && min(which(alf != 0)) < 59) {
                    FastqQuality
                } else SFastqQuality
            }, SFastqQuality=SFastqQuality, FastqQuality=FastqQuality)
        quality <- qualityFunc(quality)
        srq <- 
            if (withIds)
                ShortReadQ(sread, quality, id)
            else
                ShortReadQ(sread, quality)
        if (!missing(filter))
            srq <- srq[filter(srq)]
        srq
    }, error=function(err) {
        .throw(SRError("IncompatibleTypes", "message: %s",
                       conditionMessage(err)))
    })
})

setMethod(ShortReadQ, c("DNAStringSet", "BStringSet", "missing"),
          function(sread, quality, id, ...)
{
    ShortReadQ(sread, quality, BStringSet(character(length(sread))),
               ...)
})

setMethod(ShortReadQ, c("missing", "missing", "missing"),
          function(sread, quality, id, ...)
{
    new("ShortReadQ")
})

setAs("ShortReadQ", "QualityScaledDNAStringSet", function(from) 
{
    q <- quality(from)
    q <-
        if (is(q, "SFastqQuality")) as(q, "SolexaQuality")
        else if (is(q, "FastqQuality")) as(q, "PhredQuality")
        else as(q, "XStringQuality")
    QualityScaledDNAStringSet(sread(from), q)
})

setMethod(readFastq, "character",
    function(dirPath, pattern=character(0), ...,
             withIds=TRUE)
{
    src <- .file_names(dirPath, pattern)
    tryCatch({
        elts <- .Call(.read_solexa_fastq, src, withIds)
        if (withIds)
            ShortReadQ(elts[["sread"]], elts[["quality"]], elts[["id"]], ...,
                       withIds=withIds)
        else
            ShortReadQ(elts[["sread"]], elts[["quality"]], ...,
                       withIds=withIds)
    }, error=function(err) {
        .throw(SRError("Input/Output",
                       "file(s):\n    %s\n  message: %s",
                       paste(src, collapse="\n    "),
                       conditionMessage(err)))
    })
})

setMethod(writeFastq, c("ShortReadQ", "character"),
    function(object, file, mode="w", full=FALSE, ...)
{
    if (length(file) != 1)
        .throw(SRError("UserArgumentMismatch", "'%s' must be '%s'",
                       "file", "character(1)"))
    if (file.exists(file) && mode != "a")
        .throw(SRError("UserArgumentMismatch",
                       "file '%s' exists, but mode is not 'a'",
                       file))
    file <- path.expand(file)
    ## FIXME: different quality types
    max_width <- max(c(unique(width(id(object))),
                       unique(width(sread(object))),
                       unique(width(quality(object)))))
    if (!is(quality(quality(object)), "XStringSet"))
        .throw(SRError("UserArgumentMismatch", "'is(<%s>, \"%s\")' failed",
                       "quality", "XStringSet"))
    .Call(.write_fastq, id(object), sread(object),
          quality(quality(object)), file, mode, full, max_width)
    invisible(length(object))
})

## coerce

setMethod(pairwiseAlignment, "ShortReadQ",
          function(pattern, subject, ...)
          {
            mc <- as.list(match.call())
            if (is.null(mc$patternQuality))
              mc$patternQuality <- quality(quality(pattern))
            do.call(callNextMethod, c(list(pattern, subject), mc))
          })

## subset

setMethod("[", c("ShortReadQ", "missing", "missing"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("ShortReadQ", "missing", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("ShortReadQ", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

.ShortReadQ_subset <- function(x, i, j, ..., drop=TRUE) {
    if (0L != length(list(...))) .subset_err()
    initialize(x, sread=sread(x)[i], id=id(x)[i],
               quality=quality(x)[i])
}

setMethod("[", c("ShortReadQ", "ANY", "missing"), .ShortReadQ_subset)

setMethod(append, c("ShortReadQ", "ShortReadQ", "missing"),
    function(x, values, after=length(x))
{
    initialize(x, id=append(id(x), id(values)),
               sread=append(sread(x), sread(values)),
               quality=append(quality(x), quality(values)))
})

setMethod(narrow, "ShortReadQ",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
{
    initialize(x,
               sread=narrow(sread(x), start, end, width, use.names),
               quality=narrow(quality(x), start, end, width, use.names))
})

## manip

.abc_ShortReadQ <- function(stringSet, alphabet, ...)
{
    if (!missing(alphabet)) {
        if (!(is.list(alphabet) && length(alphabet) == 2))
            .throw(SRError("UserArgumentMismatch",
                           "'%s' must be '%s'", "alphabet",
                           "list(2)"))
        if (!all(sapply(alphabet, is, "character")))
            .throw(SRError("UserArgumentMismatch",
                           "'%s' list elements must be '%s'",
                           "alphabet", "character()"))
    }
    sread <- sread(stringSet)
    quality <- quality(stringSet)
    if (missing(alphabet))
        alphabet <- list(Biostrings::alphabet(sread),
                         Biostrings::alphabet(quality))
    w <- max(0L, width(stringSet))
    res <- .Call(.alphabet_pair_by_cycle, sread, quality(quality),
                 w, alphabet[[1]], alphabet[[2]])
    dm <- dimnames(res)
    dm[[3]]<- seq_len(w)
    names(dm)[[3]] <- "cycle"
    dimnames(res) <- dm
    res
}

setMethod(alphabetByCycle, "ShortReadQ", .abc_ShortReadQ)

setMethod(alphabetScore, "ShortReadQ", .forward_objq)

setMethod(trimTailw, "ShortReadQ",
    function(object, k, a, halfwidth, ..., ranges=FALSE)
{

    rng <- callGeneric(quality(object), k, a, halfwidth, ...,
                       ranges=TRUE)
    if (ranges) rng
    else narrow(object, 1L, end(rng))[0L != width(rng)]
})

setMethod(trimTails, "ShortReadQ",
    function(object, k, a, successive=FALSE, ..., ranges=FALSE)
{
    rng <- callGeneric(quality(object), k, a, successive,
                       ..., ranges=TRUE)
    if (ranges) rng
    else narrow(object, 1L, end(rng))[0L != width(rng)]
})

setMethod(trimEnds, "ShortReadQ",
    function(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
             ..., ranges=FALSE)
{
    rng <- callGeneric(quality(object), a, left, right, relation,
                       ..., ranges=TRUE)
    if (ranges) rng
    else narrow(object, start(rng), end(rng))
})

## show

setMethod(detail, "ShortReadQ", function(x, ...) {
    callNextMethod()
    detail(quality(x))
})

## summary

## perhaps summary stats like ShortRead except with qualities
