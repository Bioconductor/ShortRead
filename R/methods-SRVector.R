setMethod(.srValidity, "SRVector", function(object) {
    msg <- NULL
    cls <- vclass(object)
    if (length(cls)!=1)
        msg <- c(msg, "'vclass' must be character(1)")
    if (!all(sapply(object, is, cls)))
        msg <- c(msg,
                 sprintf("all elements must satisfy 'is(element, \"%s\")'",
                         cls))
    if (is.null(msg)) TRUE else msg
})

SRVector <- function(..., vclass) {
    args <- list(...)
    if (length(args)>0 && missing(vclass))
        vclass <- class(args[[1]])
    ok <- sapply(args, is, vclass)
    if (!all(ok)) {
        classes <- paste(unique(c(sapply(args, class), vclass)),
                         collapse="' '")
        .throw(SRError("SRVectorClassDisagreement",
                      "elements and vclass: '%s'", classes),
              call=match.call())
    }
    new("SRVector", .srlist=args, vclass=vclass)
}

.make_getter("vclass")

setMethod(show, "SRVector", function(object) {
    callNextMethod()
    cat("vclass: ", vclass(object), "\n", sep="")
})

setMethod(detail, "SRVector", function(x) {
    .SRList_show_class(x)
    show(unlist(.srlist(x)))
})
