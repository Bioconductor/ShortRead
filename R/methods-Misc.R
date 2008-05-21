.abc_BStringSet <- function(stringSet, alphabet, ...) {
    if (missing(alphabet))
        alphabet <- sapply(as.raw(1:255), rawToChar)
    callNextMethod(stringSet, alphabet=alphabet)
}

setMethod("alphabetByCycle", "BStringSet", .abc_BStringSet)

setMethod("srorder", "XStringSet", function(x, ...) {
    if (length(list(...))!=0)
        .throw(SRError("UserArgumentMismatch",
                       "argument '%s' not supported",
                       names(list(...))))
    .Call(.alphabet_order, x)
})

setMethod("srsort", "XStringSet", function(x, ...) x[srorder(x, ...)])

.alf_srduplicated <- function(x, incomparables=FALSE, ...) {
    if (!missing(incomparables))
        .throw(SRError("UserArgumentMismatch",
                       "argument '%s' not supported", "incomparables"))
    .Call(.alphabet_duplicated, x)
}

setMethod("srduplicated", "XStringSet", .alf_srduplicated)
