setMethod(.srValidity, "SolexaSet", function(object) {
    msg <- NULL
    nr <- nrow(laneDescription(object))
    if (nr!=8)
        msg <- c(msg,
                 sprintf("'laneDescription' must have 8 rows, but has %d",
                         nr))
    if (is.null(msg)) TRUE else msg
})

.SolexaSet_SolexaPath <- function(path, laneDescription, ...) {
    if (missing(laneDescription)) {
        laneDescription <-
            new("AnnotatedDataFrame", data=data.frame(1:8)[,FALSE],
                varMetadata=data.frame(labelDescription=character(0)),
                dimLabels=c("laneNames", "laneColumns"))
    } else {
        if (!is(laneDescription, "AnnotatedDataFrame")) {
            cls <- paste(class(laneDescription), collapse=" ")
            .throw(SRError("UserArgumentMismatch",
                          "expected '%s' as '%s', but got '%s'",
                          "AnnotatedDataFrame",
                          "laneDescription", cls))
        }
        dimLabels(laneDescription) <- c("laneNames", "laneColumns")
    }

    new("SolexaSet", ..., solexaPath=path,
        laneDescription=laneDescription)
}

setMethod(SolexaSet, "SolexaPath", .SolexaSet_SolexaPath)

setMethod(SolexaSet, "character", function(path, ...) {
    .SolexaSet_SolexaPath(SolexaPath(path), ...)
})

.make_getter(slotNames("SolexaSet"))

setMethod(laneNames, "SolexaSet", function(object, ...) {
    laneNames(laneDescription(object))
})

setMethod(laneNames, "AnnotatedDataFrame", function(object) {
    sampleNames(object)
})

## .qa_SolexaSet <- function(dirPath, pattern=character(0), ...)
## {
##     dirPath <- analysisPath(dirPath)
##     if (missing(pattern))
##         pattern <- ".*_export.txt$"
##     callGeneric(dirPath, pattern, type="SolexaExport", ...)
## }

## setMethod(qa, "SolexaSet", .qa_solexa_export)

## .report_SolexaSet <- function(x, run=1, ..., qaFile=tempfile(),
##                               dest=tempfile(), type="pdf" )
## {
##     report(qa(x, run=run))
## }

## alignment

.readAligned_SolexaSet <- function(dirPath,
                                   pattern=".*_export.txt$",
                                   run, ...)
{
    dirPath <- analysisPath(solexaPath(dirPath))[run]
    .readAligned_character(dirPath, pattern, type="SolexaExport", ...)
}

setMethod(readAligned, "SolexaSet", .readAligned_SolexaSet)

setMethod(show, "SolexaSet", function(object) {
    callNextMethod()
    cat("experimentPath(solexaPath(object)):\n  ",
        experimentPath(solexaPath(object)), "\n", sep="")
    cat("laneDescription:\n")
    print(laneDescription(object))
})

setMethod(detail, "SolexaSet", function(object, ...) {
    callNextMethod()
    cat("\n")
    detail(solexaPath(object), ...)
    cat("\nclass: AnnotatedDataFrame\n")
    ld <- laneDescription(object)
    cat("pData:\n")
    print(pData(ld))
    cat("varMetadata:\n")
    print(varMetadata(ld))
})
