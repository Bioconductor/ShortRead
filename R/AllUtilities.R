## Errors

.undefined_method_err <- function(class, method) {
  .throw(SRError("InternalError",
                 "undefined method '%s' for class '%s'",
                 method, class))
}

.subset_err <- function() {
    .throw(SRError("UserSubset",
                   "'[' must be called with only subscript 'i'"))
}

.arg_mismatch_type_err <- function(arg, type) {
    .throw(SRError("UserArgumentMismatch", "'%s' must be '%s'",
                   arg, type))
}

## Misc

.order_chr_levels <- function(lvls) {
    idx <- grep("chr[0-9XYM]+$", lvls)
    num <- grep("chr[0-9]+$", lvls[idx])
    numIds <- lvls[idx][num]
    c(numIds[order(as.numeric(sub("chr", "", numIds)))],
      lvls[idx][-num], lvls[-idx])
}

.file_names <- function(dirPath, pattern, ..., full.names=TRUE) {
    if (!is(pattern, "character") || length(pattern)>1)
        .arg_mismatch_type_err("pattern", "character(0) or character(1)")
    if (!isTRUE(full.names))
        .arg_mismatch_type_err("full.names", "TRUE")
    files <- list.files(path.expand(dirPath), pattern,
                        ..., full.names=full.names)
    files <- files[!file.info(files)$isdir]
    if (length(files)==0) {
        if (length(pattern)==0) pattern <- "character(0)"
        .throw(SRError("Input/Output",
                       "no input files found\n  dirPath: %s\n  pattern: %s\n",
                       dirPath, pattern))
    }
    files
}

.show_some <- function(what, obj) {
    if (length(obj) == 0)
      cat(what, ": (0 total)\n", sep="")
    else
      cat(what, ": ", paste(selectSome(obj), collapse=" "),
          " (", length(obj), " total)\n", sep="")
}

## Class- and method-related

.forward_objq <- function(object, ...)
    callGeneric(quality(object), ...)

.forward_xq <- function(x, ...) callGeneric(quality(x), ...)

.forward_obj <- function(object, ...)
    callGeneric(sread(object), ...)

.forward_x <- function(x, ...) callGeneric(sread(x), ...)

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
