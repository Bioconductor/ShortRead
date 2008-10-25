## public

polyn <- function(nucleotides, n)
{
    if (!is.character(nucleotides) || length(nucleotides)==0)
        .throw(SRError("UserArgumentMismatch", "'%s' must be '%s'",
                       "nucleotides", "character(1) or longer"))
    if (!all(sapply(nucleotides, nchar) == 1))
        .throw(SRError("UserArgumentMismatch",
                       "'%s' must all have %d charactecters",
                       "nucleotides", 1))
    if (!is.numeric(n) || length(n) != 1)
        .throw(SRError("UserArgumentMismatch", "'%s' must be '%s'",
                       "n", "numeric(1)"))
    sapply(nucleotides,
           function(elt) paste(rep(elt, n), collapse=""))
}

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

.arg_missing_err <- function(arg, method, help) {
    .throw(SRError("UserArgumentMismatch",
                   "argument '%s' required for '%s'\n  see %s",
                   arg, method, help))
}

.arg_mismatch_type_err <- function(arg, type) {
    .throw(SRError("UserArgumentMismatch",
                   "'%s' must be '%s'",
                   arg, type))
}

.arg_mismatch_type_err2 <- function(arg, type, was) {
    .throw(SRError("UserArgumentMismatch",
                   "'%s' must be '%s', was '%s'",
                   arg, type, was))
}

.arg_mismatch_value_err <- function(arg, value, possible_vals) {
    .throw(SRError("UserArgumentMismatch",
                   "arugment '%s' had value '%s'\n  allowable values: '%s'",
                   arg, value,
                   paste(possible_vals, collapse="' '")))
}

.check_type_and_length <- function(x, type, len)
{
    name <- deparse(substitute(x))
    if (!is(x, type))
        .arg_mismatch_type_err2(name, type, class(x))
    if (!is.na(len) && sum(length(x) == len)==0) {
        typelen <- paste(type, paste("(", len, ")", sep=""),
                         sep="", collapse="' '")
        was <- sprintf("%s(%d)", class(x), length(x))
        .arg_mismatch_type_err2(name, typelen, was)
    }
}

## Misc

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

.make_getter <-
    function(slots, where=topenv(parent.frame()), verbose=FALSE)
{
    slots <- .nameAll(slots)
    nms <- names(slots)
    ok <- !sapply(nms, exists, where)
    if (verbose && !all(ok))
        .throw(SRError("InternalError",
                      "getter '%s' already exists",
                      paste(nms[!ok], collapse=", ")))
    slots <- slots[ok]
    for (i in seq_along(slots)) {
        func <- eval(substitute(function(object, ...) slot(object, SLOT),
                                list(SLOT=slots[i])))
        assign(nms[i], func, where)
    }
}
