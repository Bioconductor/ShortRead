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

.qa_lst_as_data_frame <- function(lst) {
    if (length(lst)==0) return(data.frame())
    nms <- names(lst[[1]])
    sublst <- sapply(nms, function(nm) {
        subListExtract(lst, nm, simplify=TRUE)
    })
    names(sublst) <- nms
    do.call("data.frame", sublst)
}

.qa_Solexa_tileStats <- function(dirPath, pattern, ...) {
    .qa_Solexa_tileStats_tile <- function(dirPath, pattern, ...) {
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
    lst <- srapply(list.files(dirPath, pattern),
                   .qa_Solexa_tileStats_tile, dirPath=dirPath)
    .qa_lst_as_data_frame(lst)
}

.qa_solexa_export <- function(dirPath, pattern, type="SolexaExport", ...) {
    .lane <- function(dirPath, pattern, ..., verbose=FALSE) {
        rpt <- readAligned(dirPath, pattern, ...)
        alf <- alphabetFrequency(sread(rpt), baseOnly=TRUE)
        cleanIdx <- alf[,'other']==0
        df <- pData(alignData(rpt))
        filterIdx <- df$filtering=="Y"
        mapIdx <- !is.na(position(rpt))

        rm(alf)

        nReadByTile <- table(df$tile)
        nFilterByTile <- table(df$tile[filterIdx])
        nCleanByTile <- table(df$tile[cleanIdx])
        nFilterCleanByTile <- table(df$tile[filterIdx & cleanIdx])
        nMapByTile <- table(df$tile[mapIdx])

        baseByCycle <- alphabetByCycle(sread(rpt))
        qualityByCycle <- alphabetByCycle(quality(rpt))

        qualityScore <- alphabetScore(quality(rpt)) / width(quality(rpt))
        quality <- density(qualityScore)
        qualityFilter <- density(qualityScore[filterIdx])
        
        aquality <- quality(alignQuality(rpt))
        alignQuality <- table(aquality[filterIdx & !is.na(aquality)])

        c1 <- as.character(sread(rpt))
        tAll <- sort(table(c1), decreasing=TRUE)
        ttAll <- table(tAll)
        tAll <- head(tAll, 40)
        tFiltered <- sort(table(c1[filterIdx]), decreasing=TRUE)
        ttFiltered <- table(tFiltered)
        tFiltered=head(tFiltered, 40)
        tMapped <- sort(table(c1[mapIdx]), decreasing=TRUE)
        ttMapped <- table(tMapped)
        tMapped <- head(tMapped, 40)

        list(lane=list(
               nRead=sum(nReadByTile),
               nFilter=sum(nFilterByTile),
               nClean=sum(nCleanByTile),
               nFilterClean=sum(nFilterCleanByTile),
               nMap <- sum(nMapByTile),
               ##
               baseCall=rowSums(baseByCycle),
               qualityScore=rowSums(qualityByCycle),
               quality=quality,
               qualityFilter=qualityFilter,
               alignQuality=alignQuality,
               readFreq=list(
                 all=tAll,
                 tAll=ttAll,
                 filter=tFilter,
                 tFilter=ttFilter,
                 map=tMap,
                 tMap=ttMap),
               perCycle=list(
                 baseCall=baseByCycle,
                 qualityScore=qualityByCycle)),
             perTile=list(
               nRead=nReadByTile,
               nFilter=nFilterByTile,
               nClean=nCleanByTile,
               nFilterClean=nFilterCleanByTile,
               nMap <- nMapByTile,

               qualityRead=tapply(qualityScore, df$tile, median),
               qualityFilter=tapply(qualityScore[filterIdx],
                 df$tile[filterIdx], median),
               qualityClean=tapply(qualityScore[cleanIdx], df$tile[cleanIdx], median),
               qualityFilterClean=tapply(qualityScore[filterIdx & cleanIdx],
                 df$tile[filterIdx & cleanIdx], median)))
    }
    fls <- list.files(dirPath, pattern, full.names=TRUE)
    lst <- srapply(basename(fls), .lane, dirPath=dirPath, type=type)
    names(lst) <- basename(fls)
    lanes <- lapply(lst, "[[", "lane")

    .densityPlot <- function(lanes, part) {
        x <- lapply(lanes, function(elt) elt[[part]]$x)
        y <- lapply(lanes, function(elt) elt[[part]]$y)
        qualityDf <- data.frame(quality=unlist(x, use.names=FALSE),
                                density=unlist(y, use.names=FALSE),
                                name=rep(names(lanes), sapply(x, length, USE.NAMES=FALSE)))
        xyplot(density~quality|name, qualityDf, type="l")
    }

    list(counts=sapply(c("nRead", "nFilter", "nClean",
           "nFilterClean", "nMap"),
           function(elt) sapply(lanes, "[[", elt)),
         baseCall=t(sapply(lanes, "[[", "baseCall")),
         qualityRead=.densityPlot(lanes, "quality"),
         qualityFilter=.densityPlot(lanes, "qualityFilter"),
         readFreq=lapply(lanes, "[[", "readFreq"))
}

.qa_character <- function(dirPath, pattern=character(0),
                          type=c("SolexaExport"), ...) {
    tryCatch(type <- match.arg(type),
             error=function(err) {
                 .throw(SRError("UserArgumentMismatch",
                                conditionMessage(err)))
             })
    switch(type,
           SolexaExport=.qa_solexa_export(dirPath, pattern,
             type="SolexaExport", ...))
}

setMethod("qa", "character", .qa_character)

setMethod("qa", "SolexaSet", .qa_solexa_export)

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
