.nameAll <- function(x) {
    ## Add names to character vector x.  Elements of x without names get
    ## a name matching the element.
    if (!is.character(x))
      stop("argument 'x' must be a character vector")
    if (length(names(x)))
      names(x) <- ifelse(nchar(names(x)) == 0, x, names(x))
    else
      names(x) <- x
    x
}

.make_getter <- function(slots, where=topenv(parent.frame())) {
    slots <- .nameAll(slots)
    nms <- names(slots)
    ok <- !sapply(nms, exists, where)
    if (!all(ok))
        .throw(SRError("InternalError",
                      "'getter' already exists: %s",
                      paste(nms[!ok], collapse=", ")))
    for (i in seq_along(slots)) {
        func <- eval(substitute(function(object, ...) slot(object, SLOT),
                                list(SLOT=slots[i])))
        assign(nms[i], func, where)
    }
}
