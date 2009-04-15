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

setMethod(ShortReadQ, c("missing", "missing", "missing"),
          function(sread, quality, id, ...)
{
    new("ShortReadQ")
})

setMethod(readFastq, "character",
    function(dirPath, pattern=character(), ..., withIds=TRUE) 
{
    src <- .file_names(dirPath, pattern)
    elts <- .Call(.read_solexa_fastq, src, withIds)
    quality <- SFastqQuality(elts[["quality"]])
    if (withIds)
        ShortReadQ(elts[["sread"]], quality, elts[["id"]])
    else
        ShortReadQ(elts[["sread"]], quality)
})

setMethod(writeFastq, "ShortReadQ", function(object, file, mode="w", ...) {
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
    .Call(.write_fastq, id(object), sread(object),
          quality(quality(object)), file, mode, max_width)
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
    if (nargs() != 2) .subset_err()
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

.abc_ShortReadQ <- function(stringSet, alphabet, ...) {
    if (length(stringSet)==0)
        .throw(SRError("UserArgumentMismatch",
                       "'stringSet' must have non-zero length"))
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
    res <- .Call(.alphabet_pair_by_cycle, sread, quality(quality),
                 unique(width(stringSet)), alphabet[[1]], alphabet[[2]])
    dm <- dimnames(res)
    dm[[3]]<- seq_len(unique(width(stringSet)))
    names(dm)[[3]] <- "cycle"
    dimnames(res) <- dm
    res
}

setMethod(alphabetByCycle, "ShortReadQ", .abc_ShortReadQ)

setMethod(alphabetScore, "ShortReadQ", .forward_objq)

## show

setMethod(detail, "ShortReadQ", function(object, ...) {
    callNextMethod()
    detail(quality(object))
})

## summary

## perhaps summary stats like ShortRead except with qualities
