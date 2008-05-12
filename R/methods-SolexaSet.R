setMethod(".srValidity", "SolexaSet", function(object) {
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

setMethod("SolexaSet", "SolexaPath", .SolexaSet_SolexaPath)

setMethod("SolexaSet", "character", function(path, ...) {
    .SolexaSet_SolexaPath(SolexaPath(path), ...)
})

.make_getter(slotNames("SolexaSet"))

setMethod("laneNames", "SolexaSet", function(object, ...) {
    laneNames(laneDescription(object))
})

setMethod("laneNames", "AnnotatedDataFrame", function(object) {
    sampleNames(object)
})

## qa

.qa_SolexaSet_readCount <- function(bcPath, pattern, ...) {
    .qa_SolexaSet_readCount_tiles <- function(dirPath, pattern, ...) {
        lane <- as.numeric(sub("s_([0-9]+)_.*", "\\1", pattern))
        tile <- as.numeric(sub("s_[0-9]+_([0-9]+)_.*", "\\1", pattern))
        dna <- readXStringColumns(dirPath, pattern,
                                  colClasses=list(NULL, NULL, NULL, NULL,
                                    "DNAString"))[[1]]
        list(lane=lane, tile=tile,
             slane=(lane-1)*3+trunc((tile-1)/100)+1,
             stile=1+pmin((tile-1)%%200, (200-tile)%%200),
             nReads=length(dna),
             nClean=sum(alphabetFrequency(dna, baseOnly=TRUE)[,"other"]==0))
    }
    if (length(pattern)==0) pattern=".*_seq.txt"
    res <- srapply(list.files(bcPath, pattern),
                   .qa_SolexaSet_readCount_tiles, dirPath=bcPath)
    sublst <- lapply(names(res[[1]]), function(nm) {
        subListExtract(res, nm, simplify=TRUE)
    })
    names(sublst) <- names(res[[1]])
    do.call("data.frame", sublst)
}

.qa_SolexaSet <- function(set, pattern=character(0), baseCallRun=1, ...) {
    bcPath <- baseCallPath(solexaPath(set))[[baseCallRun]]
    list(readCount=.qa_SolexaSet_readCount(bcPath, pattern=pattern, ...))
}

setMethod("qa", "SolexaSet", .qa_SolexaSet)

## alignment

.readAligned_SolexaSet <- function(dirPath,
                                   pattern="s_[1-8]_export.txt",
                                   run=1, ...) {
    dirPath <- analysisPath(solexaPath(dirPath))[[run]]
    .readAligned_SolexaExport(dirPath, pattern, ...)
}

setMethod("readAligned", "SolexaSet", .readAligned_SolexaSet)

setMethod("show", "SolexaSet", function(object) {
    callNextMethod()
    cat("experimentPath(solexaPath(object)):\n  ",
        experimentPath(solexaPath(object)), "\n", sep="")
    cat("laneDescription:\n")
    print(laneDescription(object))
})

setMethod("detail", "SolexaSet", function(object, ...) {
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
