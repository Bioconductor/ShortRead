.abc_BStringSet <- function(stringSet, alphabet, ...) {
    if (missing(alphabet))
        alphabet <- sapply(as.raw(0:255), rawToChar)
    callNextMethod(stringSet, alphabet=alphabet)
}

setMethod("alphabetByCycle", "BStringSet", .abc_BStringSet)
