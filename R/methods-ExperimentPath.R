setMethod(".srValidity", "ExperimentPath", function(object) {
    msg <- NULL
    if (length(basePath(object))!=1)
        msg <- c(msg, "ExperimentPath 'basePath' must be character(1)")
    if (is.null(msg)) TRUE else msg
})

.srPath <- function(path, pattern) {
    path <- path.expand(path)
    tryCatch({
        res <- list.files(path, pattern=pattern, full.name=TRUE)
        if (length(res)==0) NA_character_
        else res
    }, warning=function(warn) NA_character_)
}

.checkPath <- function(path) {
  nm <- deparse(substitute(path))
  if (length(path)==0) {
    warning(nm, " not defined")
  } else {
    for (p in path)
      if (!file.exists(p)) 
        warning(nm, " '", p, "' does not exist")
  }
}

.make_getter(slotNames("ExperimentPath"))

setMethod("show", "ExperimentPath", function(object) {
    catPath <- function(nm) {
        vals <- do.call(nm, list(object))
        vals <- substr(basename(vals), 1, 15)
        vals <- paste(vals, ifelse(nchar(vals)==15, "...", ""),
                      sep="")
        cat(nm, ": ", paste(vals, collapse=", "), "\n", sep="")
    }
    callNextMethod()
    cat("basePath: ", basePath(object), "\n", sep="")
    slts <- slotNames(object)
    for (slt in slts[slts!="basePath"]) catPath(slt)
})

setMethod("detail", "ExperimentPath", function(object, ...) {
    catPath <- function(nm) {
        fnms <- do.call(nm, list(object))
        cat(nm, ":\n  ", paste(fnms, collapse="\n  "), sep="")
        cat("\n")
    }
    callNextMethod()
    cat("basePath:\n  ", basePath(object), "\n", sep="")
    slts <- slotNames(object)
    for (slt in slts[slts!="basePath"]) catPath(slt)
})
