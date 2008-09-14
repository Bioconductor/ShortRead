setMethod(".srValidity", "SolexaPath", function(object) {
    msg <- NULL
    if (length(experimentPath(object))!=1)
        msg <- c(msg, "SolexaPath 'experimentPath' must be character(1)")
    if (is.null(msg)) TRUE else msg
})

SolexaPath <- function(experimentPath=NA_character_,
                       dataPath=.srPath(experimentPath, "Data"),
                       scanPath=.srPath(dataPath, "GoldCrest"),
                       imageAnalysisPath=.srPath(dataPath, "^C"),
                       baseCallPath=.srPath(imageAnalysisPath,
                         "^Bustard"),
                       analysisPath=.srPath(baseCallPath,
                         "^GERALD"),
                       ..., verbose=FALSE) {
    if (verbose) {
        .checkPath(experimentPath)
        .checkPath(dataPath)
        .checkPath(scanPath)
        .checkPath(imageAnalysisPath)
        .checkPath(baseCallPath)
        .checkPath(analysisPath)
    }
    new("SolexaPath", ..., basePath=experimentPath, dataPath=dataPath,
        scanPath=scanPath, imageAnalysisPath=imageAnalysisPath,
        baseCallPath=baseCallPath, analysisPath=analysisPath)
}

.make_getter(slotNames("SolexaPath"))

.readPrb_SolexaPath <- function(dirPath, pattern, run, ...)
{
    callGeneric(baseCallPath(dirPath)[run], pattern, ...)
}

setMethod("readPrb", "SolexaPath", .readPrb_SolexaPath)

.readFastq_SolexaPath <- function(dirPath, 
                                  pattern=".*_sequence.txt$",
                                  run, ...) {
    dirPath <- analysisPath(dirPath)[run]
    if (is.na(dirPath))
        .throw(SRError("Input/Output", "'%s' is 'NA' in '%s'",
                      "analysisPath", "dirPath"))
    callGeneric(dirPath, ..., pattern=pattern)
}

setMethod("readFastq", "SolexaPath", .readFastq_SolexaPath)

.readBaseQuality_SolexaPath <- function(dirPath,
                                        seqPattern=".*_seq.txt$",
                                        prbPattern=".*_prb.txt$",
                                        run, ...)
{
    dirPath <- baseCallPath(dirPath)[run]
    .readBaseQuality_Solexa(dirPath, seqPattern=seqPattern,
                            prbPattern=prbPattern, ...)
}

setMethod("readBaseQuality", "SolexaPath", .readBaseQuality_SolexaPath)

.readAligned_SolexaPath <- function(dirPath,
                                    pattern=".*_export.txt$",
                                    run, ...) 
{
    dirPath <- analysisPath(dirPath)[run]
    .readAligned_SolexaExport(dirPath, pattern, ...)
}

setMethod("readAligned", "SolexaPath", .readAligned_SolexaPath)

.qa_SolexaPath <- function(dirPath, pattern=character(0),
                           run, ...)
{
    dirPath <- analysisPath(dirPath)[run]
    if (missing(pattern))
        pattern <- ".*_export.txt$"
    callGeneric(dirPath, pattern, type="SolexaExport", ...)
}

setMethod("qa", "SolexaPath", .qa_SolexaPath)

.report_SolexaPath <- function(x, ..., dest=tempfile(), type="pdf" )
{
    report(qa(x, ...), dest=dest, type=type)
}

setMethod("report", "SolexaPath", .report_SolexaPath)

setMethod("show", "SolexaPath", function(object) {
    callNextMethod()
    .show_additionalPathSlots(object)
})

setMethod("detail", "SolexaPath", function(object, ...) {
    callNextMethod()
    .detail_additionalPathSlots(object)
})
