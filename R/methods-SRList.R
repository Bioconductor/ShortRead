SRList <- function(...) {
    args <- list(...)
    if (length(args)==1 && is(args[[1]], "list"))
        new("SRList", .srlist=args[[1]])
    else
        new("SRList", .srlist=args)
}

.srlist <- .make_getter(".srlist")

setMethod("names", "SRList", function(x) names(.srlist(x)))

setReplaceMethod("names", c("SRList", "character"), function(x, value) {
    lst <- .srlist(x)
    names(lst) <- value
    initialize(x, .srlist=lst)
})

setMethod("length", "SRList", function(x) length(.srlist(x)))

setMethod("[", c(x="SRList", i="ANY", j="missing"),
          function(x, i, j, ..., drop=FALSE) {
              initialize(x, .srlist=.srlist(x)[i])
          })

setMethod("[[", signature(x="SRList", i="ANY", j="missing"),
          function(x, i, j, ...) .srlist(x)[[i]])

setMethod("sapply", "SRList", function(X, FUN, ..., simplify=TRUE,
                                       USE.NAMES=TRUE) {
    sapply(.srlist(X), FUN, ..., simplify=simplify,
           USE.NAMES=USE.NAMES)
})

setMethod("lapply", "SRList", function(X, FUN, ...) {
    lapply(.srlist(X), FUN, ...)
})

.SRList_show_class <- function(object) {
    cat("class: ", class(object), "(", length(object), ")\n", sep="")
}

setMethod("show", "SRList", .SRList_show_class)

setMethod("detail", "SRList", function(object,...) {
    .SRList_show_class(object)
    .srlist(object)
})
