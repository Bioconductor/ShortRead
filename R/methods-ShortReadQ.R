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

.make_getter("quality")

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
    initialize(x, sread=sread(x)[i], id=id(x)[i], quality=quality(x)[i])
}

setMethod("[", c("ShortReadQ", "ANY", "missing"), .ShortReadQ_subset)

## manip

.abc_ShortReadQ <- function(stringSet, alphabet, ...) {
    if (!missing(alphabet))
        .throw(SRWarn("UserArgumentMismatch", "'alphabet' ignored"))
    res <- lapply(c(sread(stringSet), quality(stringSet)),
                  alphabetByCycle, ...)
    names(res) <- c("sread", "quality")
    res
}

setMethod("alphabetByCycle", "ShortReadQ", .abc_ShortReadQ)

.ascore_ShortReadQ <- function(object, ...) {
    callGeneric(quality(object), ...)
}

setMethod("alphabetScore", "ShortReadQ", .ascore_ShortReadQ)

## show

setMethod("detail", "ShortReadQ", function(object, ...) {
    callNextMethod()
    detail(quality(object))
})
