## SRUtil

setGeneric(".throw",
           function(object, call=NULL, ...) standardGeneric(".throw"),
           signature=c("object"))

countLines <- function(dirPath, ..., pattern=character(0), useFullName=FALSE)
{
    src <- list.files(dirPath, ..., pattern=pattern, full.name=TRUE)
    if (length(src)==0)
        .throw(SRError("Input/Output",
                      "no files in directory '%s' matching '%s'",
                      dirPath, ifelse(length(pattern)==0, "", pattern)))
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

## ShortRead / ShortReadQ

setGeneric("readFastq", function(dirPath, pattern=character(0), ...)
           standardGeneric("readFastq"), signature="dirPath")

setGeneric("length")

setGeneric("clean", function(object, ...) standardGeneric("clean"))

detail <- function(object, ...) show(object)

setGeneric("detail")

## AlignedRead

setGeneric("readAligned", function(dirPath, pattern=character(0), ...)
           standardGeneric("readAligned"), signature="dirPath")

## SolexaSet

setGeneric("SolexaSet", function(path, ...)
           standardGeneric("SolexaSet"))

setGeneric("laneNames", function(object, ...) {
    standardGeneric("laneNames")
})
