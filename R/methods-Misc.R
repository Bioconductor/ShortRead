.abc_BStringSet <- function(stringSet, alphabet, ...) {
    if (missing(alphabet))
        alphabet <- sapply(as.raw(1:255), rawToChar)
    callNextMethod(stringSet, alphabet=alphabet)
}

setMethod("clean", "DNAStringSet", function(object, ...) {
    object[alphabetFrequency(object, baseOnly=TRUE)[,'other']==0]
})

setMethod("alphabetByCycle", "BStringSet", .abc_BStringSet)

setMethod("srorder", "XStringSet", function(x, ...) {
    if (length(list(...))!=0)
        .throw(SRError("UserArgumentMismatch",
                       "argument '%s' not supported",
                       names(list(...))))
    .Call(.alphabet_order, x)
})

setMethod("srrank", "XStringSet", function(x, ...) {
    if (length(list(...))!=0)
        .throw(SRError("UserArgumentMismatch",
                       "argument '%s' not supported",
                       names(list(...))))
    .Call(.alphabet_rank, x)
})

setMethod("srsort", "XStringSet", function(x, ...) x[srorder(x, ...)])

setMethod("srduplicated", "XStringSet", function(x, ...) {
    if (length(list(...))!=0)
        .throw(SRError("UserArgumentMismatch",
                       "argument '%s' not supported",
                       names(list(...))))
    .Call(.alphabet_duplicated, x)
})

.stringset_tables <- function(x, n=50, ...) {
    ## FIXME: two sorts
    srt <- srsort(x)
    r <- srrank(x)
    t <- tabulate(r)
    o <- order(t, decreasing=TRUE)
    ## n most common sequences
    top <- head(t[o], n)
    names(top) <- as.character(head(srt[o], n))
    ## overall frequency -- equivalent of table(table(sread))
    tt <- tabulate(t)
    nOccurrences <- seq_along(tt)[tt!=0]
    nReads <- tt[tt!=0]
    ## results
    list(top=top,
         distribution=data.frame(
           nOccurrences=nOccurrences,
           nReads=nReads))
}

setMethod("tables", "XStringSet", .stringset_tables)
