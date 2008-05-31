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
