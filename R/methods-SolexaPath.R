setMethod(".srValidity", "SolexaPath", function(object) {
    msg <- NULL
    if (length(experimentPath(object))!=1)
        msg <- c(msg, "SolexaPath 'experimentPath' must be character(1)")
    if (is.null(msg)) TRUE else msg
})

.solexaPath <- function(path, pattern) {
    path <- path.expand(path)
    tryCatch({
        res <- list.files(path, pattern=pattern, full.name=TRUE)
        if (length(res)==0) NA_character_
        else res
    }, warning=function(warn) NA_character_)
}

SolexaPath <- function(experimentPath,
                        dataPath=.solexaPath(experimentPath, "Data"),
                        scanPath=.solexaPath(dataPath, "GoldCrest"),
                        imageAnalysisPath=.solexaPath(dataPath, "^C"),
                        baseCallPath=.solexaPath(imageAnalysisPath,
                          "^Bustard"),
                        analysisPath=.solexaPath(baseCallPath,
                          "^GERALD"),
                       ..., verbose=FALSE) {
    checkPath <- function(path) {
        nm <- deparse(substitute(path))
        if (length(path)==0) {
            warning(nm, " not defined")
        } else {
            for (p in path)
                if (!file.exists(p)) 
                    warning(nm, " '", p, "' does not exist")
        }
    }
    if (verbose) {
        checkPath(experimentPath)
        checkPath(dataPath)
        checkPath(scanPath)
        checkPath(imageAnalysisPath)
        checkPath(baseCallPath)
        checkPath(analysisPath)
    }
    new("SolexaPath", ..., experimentPath=experimentPath,
        dataPath=dataPath, scanPath=scanPath,
        imageAnalysisPath=imageAnalysisPath, baseCallPath=baseCallPath,
        analysisPath=analysisPath)
}

.make_getter(slotNames("SolexaPath"))

.readFastq_SolexaPath <- function(dirPath, ...,
                                  pattern="s_[1-8]_sequence.txt") {
    dirPath <- analysisPath(dirPath)
    if (is.na(dirPath))
        .throw(SRError("Input/Output", "'%s' is 'NA' in '%s'",
                      "analysisPath", "dirPath"))
    callGeneric(dirPath, ..., pattern=pattern)
}

setMethod("readFastq", "SolexaPath", .readFastq_SolexaPath)

setMethod("show", "SolexaPath", function(object) {
    catPath <- function(nm) {
        vals <- do.call(nm, list(object))
        vals <- substr(basename(vals), 1, 15)
        vals <- paste(vals, ifelse(nchar(vals)==15, "...", ""),
                      sep="")
        cat(nm, ": ", paste(vals, collapse=", "), "\n", sep="")
    }
    callNextMethod()
    cat("experimentPath: ", experimentPath(object), "\n", sep="")
    slts <- slotNames("SolexaPath")
    for (slt in slts[slts!="experimentPath"]) catPath(slt)
})

setMethod("detail", "SolexaPath", function(object, ...) {
    catPath <- function(nm) {
        fnms <- do.call(nm, list(object))
        cat(nm, ":\n  ", paste(fnms, collapse="\n  "), sep="")
        cat("\n")
    }
    callNextMethod()
    cat("experimentPath:\n  ", experimentPath(object), "\n", sep="")
    slts <- slotNames("SolexaPath")
    for (slt in slts[slts!="experimentPath"]) catPath(slt)
})
