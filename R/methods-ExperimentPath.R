setMethod(".srValidity", "ExperimentPath", function(object) {
    msg <- NULL
    if (length(experimentPath(object))!=1)
        msg <- c(msg, "ExperimentPath 'experimentPath' must be character(1)")
    if (is.null(msg)) TRUE else msg
})

.srPath <- function(path, pattern = character()) {
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

ExperimentPath <- function(experimentPath=NA_character_, ...) {
    new("ExperimentPath", basePath=experimentPath, ...)
}

basePath <- function(object, ...) {
    .Deprecated("experimentPath")
    experimentPath(object, ...)
}

setMethod("sampleNames", "ExperimentPath", function(object) {
    character(0)
})

.show_additionalPathSlots <- function(object) { # for derived classes
    catPath <- function(nm) {
        vals <- do.call(nm, list(object))
        vals <- substr(basename(vals), 1, 15)
        vals <- paste(vals, ifelse(nchar(vals)==15, "...", ""),
                      sep="")
        cat(nm, ": ", paste(vals, collapse=", "), "\n", sep="")
    }
    slts <- slotNames(object)
    for (slt in slts[slts!="basePath"]) catPath(slt)
}

setMethod("show", "ExperimentPath", function(object) {
    callNextMethod()
    cat("experimentPath: ", experimentPath(object), "\n", sep="")
})

.detail_additionalPathSlots <- function(object) {
    catPath <- function(nm) {
        fnms <- do.call(nm, list(object))
        cat(nm, ":\n  ", paste(fnms, collapse="\n  "), sep="")
        cat("\n")
    }
    slts <- slotNames(object)
    for (slt in slts[slts!="basePath"]) catPath(slt)
}

setMethod("detail", "ExperimentPath", function(object, ...) {
    callNextMethod()
    cat("experimentPath:\n  ", experimentPath(object), "\n", sep="")
})
