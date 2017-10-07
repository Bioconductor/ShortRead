## new generics

setGeneric(".throw",
           function(object, call=NULL, ...) standardGeneric(".throw"),
           signature=c("object"))

setGeneric("renewable", function(x, ...) standardGeneric("renewable"))

setGeneric("renew", function(x, ...) standardGeneric("renew"))

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
    if (missing(alphabet))
        alphabet <- Biostrings::alphabet(stringSet)
    w <- max(0L, width(stringSet))
    .Call(.alphabet_by_cycle, stringSet, w, alphabet)
}
    
setGeneric("alphabetByCycle", signature="stringSet")

setGeneric("dustyScore", function(x, batchSize=NA, ...)
           standardGeneric("dustyScore"),
           signature="x")

setGeneric("srorder", function(x, ...) standardGeneric("srorder"))

setGeneric("srduplicated",
           function(x, ...) standardGeneric("srduplicated"))

setGeneric("srsort", function(x, ...) standardGeneric("srsort"))

setGeneric("srrank", function(x, ...) standardGeneric("srrank"))

setGeneric("tables", function(x, n=50, ...)
           standardGeneric("tables"),
           signature="x")

## Intensities

setGeneric("readIntensities",
           function(dirPath, pattern=character(0), ...)
           standardGeneric("readIntensities"),
           signature="dirPath")

## QualityScore

setGeneric("FastqQuality", function(quality, ...)
           standardGeneric("FastqQuality"))

setGeneric("SFastqQuality", function(quality, ...)
           standardGeneric("SFastqQuality"))

setGeneric("readPrb", function(dirPath, pattern=character(0), ...)
           standardGeneric("readPrb"), signature="dirPath")

## ShortRead / ShortReadQ

setGeneric("ShortRead", function(sread, id, ...)
           standardGeneric("ShortRead"))

setGeneric("sread", function(object, ...) standardGeneric("sread"))

setGeneric("writeFasta", function(object, file, mode="w", ...)
           standardGeneric("writeFasta"),
           signature=signature("object"))

setGeneric("ShortReadQ", function(sread, quality, id, ...)
           standardGeneric("ShortReadQ"))

setGeneric("readFastq", function(dirPath, pattern=character(0), ...)
           standardGeneric("readFastq"), signature="dirPath")

setGeneric("writeFastq",
           function(object, file, mode="w", full=FALSE, compress=TRUE, ...)
           standardGeneric("writeFastq"), signature=c("object", "file"))

setGeneric("readFasta",
           function(dirPath, pattern=character(0), ...,
                    nrec=-1L, skip=0L)
           standardGeneric("readFasta"),
           signature="dirPath")

setGeneric("readQual",
           function(dirPath, pattern=character(0), ...)
           standardGeneric("readQual"))

setGeneric("read454",
           function(dirPath, ...)
           standardGeneric("read454"))

setGeneric("readFastaQual",
           function(dirPath, ...)
           standardGeneric("readFastaQual"))

setGeneric("readBaseQuality",
           function(dirPath, ...)
           standardGeneric("readBaseQuality"))

setGeneric("readQseq", function(dirPath, pattern=character(0), ...,
                                as=c("ShortReadQ", "DataFrame",
                                  "XDataFrame"),
                                filtered=FALSE, verbose=FALSE)
           standardGeneric("readQseq"), signature="dirPath")

setGeneric("trimTails",
           function(object, k, a, successive=FALSE, ..., ranges=FALSE)
           standardGeneric("trimTails"), signature="object")

setGeneric("trimTailw",
           function(object, k, a, halfwidth, ..., ranges=FALSE)
           standardGeneric("trimTailw"), signature="object")

setGeneric("trimEnds",
           function(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
                    ..., ranges=FALSE)
           standardGeneric("trimEnds"),
           signature="object")

setGeneric("clean", function(object, ...) standardGeneric("clean"))

setGeneric("srdistance", function(pattern, subject, ...)
           standardGeneric("srdistance"), signature=c("pattern", "subject"))

setGeneric("alphabetScore",
           function(object, ...) standardGeneric("alphabetScore"))

## SRFilter

setGeneric("name", function(x, ...) standardGeneric("name"))

setGeneric("stats", function(x, ...) standardGeneric("stats"))

setGeneric("srFilter", function(fun, name=NA_character_, ...)
           standardGeneric("srFilter"),
           signature="fun")

## AlignedRead

setGeneric("readAligned",
           function(dirPath, pattern=character(0), ...)
           standardGeneric("readAligned"), signature="dirPath")

setGeneric("chromosome",
           function(object, ...) standardGeneric("chromosome"))

setGeneric("id",
           function(object, ...) standardGeneric("id"))

setGeneric("position",
           function(object, ...) standardGeneric("position"))

## ExperimentPath

experimentPath <- function(object, ...) {
    slot(object, "basePath")
}

setGeneric("experimentPath")

## *Set

setGeneric("qa", function(dirPath, ...) standardGeneric("qa"))

report <-
    function (x, ..., dest = tempfile(), type="html")
{
    func <- switch(type, html=report_html, pdf=.report_pdf,
                   .report_any)
    func(x, dest, type, ...)
}

setGeneric("report", signature="x")

.report_any <-
        function (x, dest, type, ...)
{
    .throw(SRError("UserArgumentMismatch",
                   "'%s, type=\"%s\"' not implemented for class '%s'",
                   "report", type, class(x)))
}

setGeneric("report_html", function(x, dest, type, ...)
           standardGeneric("report_html"),
           signature="x", useAsDefault=.report_any)

setGeneric(".report_pdf", function(x, dest, type, ...)
           standardGeneric(".report_pdf"),
           signature="x", useAsDefault=.report_any)

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

## Snapshot

setGeneric("Snapshot",
           function(files, range, ...) standardGeneric("Snapshot"))

setGeneric("SnapshotFunctionList",
           function(...) standardGeneric("SnapshotFunctionList"))

setGeneric("files", function(x, ...) standardGeneric("files"))

setGeneric("vrange", function(x, ...) standardGeneric("vrange"))

setGeneric("functions", function(x, ...) standardGeneric("functions"))

setGeneric("annTrack", function(x, ...) standardGeneric("annTrack"))

setGeneric("ignore.strand", function(x, ...) standardGeneric("ignore.strand"))

setGeneric("fac", function(x, ...) standardGeneric("fac"))

setGeneric("getTrellis", function(x, ...) standardGeneric("getTrellis"))

setGeneric("togglez", function(x, ...) standardGeneric("togglez"))

setGeneric("togglep", function(x, ...) standardGeneric("togglep"))

setGeneric("togglefun", function(x, name, ...) standardGeneric("togglefun"))

setGeneric("zoom", function(x, range, ...) standardGeneric("zoom"))

setGeneric("pan", function(x, ...) standardGeneric("pan"))

setGeneric("view", function(x, ...) standardGeneric("view"))

setGeneric("zi", function(x, ...) standardGeneric("zi"))

setGeneric("zo", function(x, ...) standardGeneric("zo"))

setGeneric("left", function(x, ...) standardGeneric("left"))

setGeneric("right", function(x, ...) standardGeneric("right"))

setGeneric("restore", function(x, ...) standardGeneric("restore"))

## ShortReadFile

setGeneric(".ShortReadFile", function(g, path, ...)
           standardGeneric(".ShortReadFile"), signature="path")

setGeneric("FastqFileList",
           function(..., class="FastqFile")
           standardGeneric("FastqFileList"),
           signature="...")
BiocGenerics:::apply_hotfix73465(getGeneric("FastqFileList"))

setGeneric("FastqStreamer",
           function(con, n, readerBlockSize=1e8, verbose=FALSE)
           standardGeneric("FastqStreamer"),
           signature=c("con", "n"))

setGeneric("FastqStreamerList",
           function(..., n, readerBlockSize=1e8, verbose=FALSE)
           standardGeneric("FastqStreamerList"),
           signature="...")
BiocGenerics:::apply_hotfix73465(getGeneric("FastqStreamerList"))

setGeneric("FastqSamplerList",
           function(..., n=1e6, readerBlockSize=1e8, verbose=FALSE,
                    ordered = FALSE)
           standardGeneric("FastqSamplerList"),
           signature="...")
BiocGenerics:::apply_hotfix73465(getGeneric("FastqSamplerList"))

setGeneric("yield", function(x, ...) standardGeneric("yield"))
