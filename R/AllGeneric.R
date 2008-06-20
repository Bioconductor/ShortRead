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
    width <- unique(Biostrings::width(stringSet))
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

## ShortRead / ShortReadQ

setGeneric("readFastq", function(dirPath, pattern=character(0), ...)
           standardGeneric("readFastq"), signature="dirPath")

setGeneric("length")

setGeneric("clean", function(object, ...) standardGeneric("clean"))

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

## SolexaSet

setGeneric("SolexaSet", function(path, ...)
           standardGeneric("SolexaSet"))

setGeneric("laneNames", function(object, ...) {
    standardGeneric("laneNames")
})
