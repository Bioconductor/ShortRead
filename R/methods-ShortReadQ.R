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

.abc_pair_by_cycle <- function(stringSet1, stringSet2, alphabet1, alphabet2) {
    if (length(stringSet1)==0)
        .throw(SRError("UserArgumentMismatch",
                       "'stringSet1' must have non-zero length"))
    if (length(stringSet2)==0)
        .throw(SRError("UserArgumentMismatch",
                       "'stringSet2' must have non-zero length"))
    if (missing(alphabet1))
        alphabet1 <- Biostrings::alphabet(stringSet1[[1]])
    if (missing(alphabet2))
        alphabet2 <- sapply(as.raw(33:132), rawToChar)
    width <- unique(c(unique(Biostrings::width(stringSet1)), unique(Biostrings::width(stringSet2))))
    if (length(width)!=1)
        .throw(SRError("UserArgumentMismatch",
                       "'width' must be unique, but is '%s'",
                       paste(width, collapse="', '")))
    .Call(.alphabet_pair_by_cycle, stringSet1, stringSet2, width, alphabet1, alphabet2)
}

.abc_ShortReadQ <- function(stringSet, alphabet, ...) {
    if (!missing(alphabet) && !(is.list(alphabet) && length(alphabet) == 2))
        .throw(SRWarn("UserArgumentMismatch", "'alphabet' ignored"))
    if (missing(alphabet))
        res <- list(sread = alphabetByCycle(sread(stringSet)),
                    quality = alphabetByCycle(quality(stringSet)),
                    both = .abc_pair_by_cycle(stringSet1 = sread(stringSet),
                                              stringSet2 = quality(quality(stringSet))))
    else
        res <-
            list(sread = alphabetByCycle(sread(stringSet), alphabet[[1]]),
                 quality = alphabetByCycle(quality(stringSet), alphabet[[2]]),
                 both = .abc_pair_by_cycle(stringSet1 = sread(stringSet),
                                           stringSet2 = quality(quality(stringSet)),
                                           alphabet1 = alphabet[[1]],
                                           alphabet2 = alphabet[[2]]))
    res
}

setMethod("alphabetByCycle", "ShortReadQ", .abc_ShortReadQ)

setMethod("alphabetScore", "ShortReadQ", .forward_objq)

## show

setMethod("detail", "ShortReadQ", function(object, ...) {
    callNextMethod()
    detail(quality(object))
})
