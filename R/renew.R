.renewable_of <-
    function(x, ...)
{
    cls <- names(getClass(x)@subclasses)
    sort(cls[!grepl("^\\.", cls)])
}

setMethod(renewable, "missing", function(x, ...) {
    ## classes that are renew-able
    .renewable_of(".ShortReadBase")
})

.renewable_query <-
    function(x)
{
    if (1L != length(x))
        .throw(SRError("UserArgumentMismatch",
                       "'%s' must be '%s'", "x", "character(1)"))
    if (!x %in% names(getClass(".ShortReadBase")@subclasses))
        .throw(SRError("UserArgumentMismatch",
                       "'%s' is not a renewable class", x))
    cls <- getClass(x)
    if (cls@virtual) {
        subcls <- .renewable_of(x)
        res <- lapply(subcls, .renewable_query)
        names(res) <- subcls
        res
    } else {
        getSlots(x)
    }
}

setMethod(renewable, ".ShortReadBase", function(x, ...) {
    structure(list(.renewable_query(class(x))), .Names=class(x))
})

setMethod(renewable, "character", function(x, ...) {
    res <- .renewable_query(x)
    if (!is.list(res))
        res <- structure(list(res), .Names=x)
    res
})

setMethod(renew, ".ShortReadBase", function(x, ...) {
    initialize(x, ...)
})
