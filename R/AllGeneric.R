## SRUtil

setGeneric(".throw",
           function(object, call=NULL, ...) standardGeneric(".throw"),
           signature=c("object"))

countLines <- function(dirPath, pattern=character(0), ..., useFullName=FALSE)
{
    src <- .file_names(path.expand(dirPath), pattern, ...)
    nLines <- .Call(.count_lines, src)
    names(nLines) <- 
        if (useFullName) src
        else basename(src)
    nLines
}

setGeneric("countLines", signature="dirPath")

alphabetByCycle <-
    function(stringSet, alphabet, ...)
{
    if (length(stringSet)==0)
        .throw(SRError("UserArgumentMismatch",
                      "'stringSet' must have non-zero length"))
    if (missing(alphabet))
        alphabet <- Biostrings::alphabet(stringSet[[1]])
    width <- unique(IRanges::width(stringSet))
    if (length(width)!=1)
        .throw(SRError("UserArgumentMismatch",
                      "'width' must be unique, but is '%s'",
                      paste(width, collapse="', '")))
    .Call(.alphabet_by_cycle, stringSet, width, alphabet)
}
    
setGeneric("alphabetByCycle", signature="stringSet")

setGeneric("srorder", function(x, ...) standardGeneric("srorder"))

setGeneric("srduplicated",
           function(x, ...) standardGeneric("srduplicated"))

setGeneric("srsort", function(x, ...) standardGeneric("srsort"))

setGeneric("srrank", function(x, ...) standardGeneric("srrank"))

setGeneric("tables", function(x, n=50, ...)
           standardGeneric("tables"),
           signature=c("x"))

## QualityScore

setGeneric("FastqQuality", function(quality, ...)
           standardGeneric("FastqQuality"))

setGeneric("SFastqQuality", function(quality, ...)
           standardGeneric("SFastqQuality"))

setGeneric("readPrb", function(dirPath, pattern=character(0), ...)
           standardGeneric("readPrb"), signature="dirPath")

## ShortRead / ShortReadQ

setGeneric("readFastq", function(dirPath, pattern=character(0), ...)
           standardGeneric("readFastq"), signature="dirPath")

setGeneric("readFasta",
           function(dirPath, pattern=character(0), ...)
           standardGeneric("readFasta"))

setGeneric("readQual",
           function(dirPath, pattern=character(0), ...)
           standardGeneric("readQual"))

setGeneric("read454",
           function(dirPath, ...)
           standardGeneric("read454"))

setGeneric("length")

setGeneric("clean", function(object, ...) standardGeneric("clean"))

setGeneric("srdistance", function(pattern, subject, ...)
           standardGeneric("srdistance"), signature=c("pattern", "subject"))

detail <- function(object, ...) show(object)

setGeneric("detail")

setGeneric("alphabetScore", function(object, ...) {
    standardGeneric("alphabetScore")
})

## AlignedRead

setGeneric("readAligned", function(dirPath, pattern=character(0), ...)
           standardGeneric("readAligned"), signature="dirPath")

## *Set

setGeneric("qa", function(dirPath, ...) standardGeneric("qa"))

setGeneric("report", function(x, ..., dest=tempfile(), type="pdf") {
    standardGeneric("report")
}, signature="x")

## SolexaSet

setGeneric("SolexaSet", function(path, ...)
           standardGeneric("SolexaSet"))

setGeneric("laneNames", function(object, ...) {
    standardGeneric("laneNames")
})

## Roche

setGeneric("RocheSet", function(path, ...)
           standardGeneric("RocheSet"))

setGeneric("runNames", function(object, ...)
           standardGeneric("runNames"))
