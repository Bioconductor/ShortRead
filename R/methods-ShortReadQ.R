## validity / accessors / constructors

setMethod(".srValidity", "ShortReadQ", function(object) {
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

setMethod("readFastq", "character", function(dirPath, pattern=character(),
                                             ...) {
    src <- .file_names(dirPath, pattern)
    elts <- .Call(.read_solexa_fastq, src)
    new("ShortReadQ", ..., sread=elts[["sread"]], id=elts[["id"]],
        quality=SFastqQuality(elts[["quality"]]))
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
                 width(stringSet), alphabet[[1]], alphabet[[2]])
    dm <- dimnames(res)
    dm[[3]]<- 1:width(stringSet)
    names(dm)[[3]] <- "cycle"
    dimnames(res) <- dm
    res
}

setMethod("alphabetByCycle", "ShortReadQ", .abc_ShortReadQ)

setMethod("alphabetScore", "ShortReadQ", .forward_objq)

## show

setMethod("detail", "ShortReadQ", function(object, ...) {
    callNextMethod()
    detail(quality(object))
})

## summary

## perhaps summary stats like ShortRead except with qualities
