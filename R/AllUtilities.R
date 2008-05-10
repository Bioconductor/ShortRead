.undefined_method_err <- function(class, method) {
  .throw(SRError("InternalError",
                 "undefined method '%s' for class '%s'",
                 method, class))
}

.subset_err <- function() {
    .throw(SRError("UserSubset",
                   "'[' must be called with only subscript 'i'"))
}

.show_some <- function(what, obj) {
    if (length(obj) == 0)
      cat(what, ": (0 total)\n", sep="")
    else
      cat(what, ": ", paste(selectSome(obj), collapse=" "),
          " (", length(obj), " total)\n", sep="")
}

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
                      "getter '%s' already exists",
                      paste(nms[!ok], collapse=", ")))
    for (i in seq_along(slots)) {
        func <- eval(substitute(function(object, ...) slot(object, SLOT),
                                list(SLOT=slots[i])))
        assign(nms[i], func, where)
    }
}
