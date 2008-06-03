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
    lst <- srapply(list.files(dirPath, pattern, full.names=TRUE),
                   .qa_Solexa_tileStats_tile, dirPath=dirPath)
    .qa_lst_as_data_frame(lst)
}

.qa_solexa_export <- function(dirPath, pattern, type="SolexaExport", ...) {
    .lane <- function(dirPath, pattern, ..., verbose=FALSE) {
        .cleanIdx <- function(rpt) {
            alphabetFrequency(sread(rpt), baseOnly=TRUE)[,'other']==0
        }
        rpt <- readAligned(dirPath, pattern,...)
        df <- pData(alignData(rpt))
        cleanIdx <- .cleanIdx(rpt)
        filterIdx <- df$filtering=="Y"
        mapIdx <- !is.na(position(rpt))

        nbins <- max(df$tile)
        nReadByTile <- tabulate(df$tile, nbins)
        nFilterByTile <- tabulate(df$tile[filterIdx], nbins)
        nCleanByTile <- tabulate(df$tile[cleanIdx], nbins)
        nFilterCleanByTile <- tabulate(df$tile[filterIdx & cleanIdx], nbins)
        nMapByTile <- tabulate(df$tile[mapIdx], nbins)

        baseByCycle <- alphabetByCycle(sread(rpt))
        qualityByCycle <- alphabetByCycle(quality(rpt))

        qualityScore <- alphabetScore(quality(rpt)) / width(quality(rpt))
        quality <- density(qualityScore)
        qualityFilter <- density(qualityScore[filterIdx])
        
        aquality <- quality(alignQuality(rpt))
        alignQuality <- table(aquality[filterIdx & !is.na(aquality)])

        dupTables <- tables(sread(rpt))
        mapDupTables <- tables(sread(rpt)[mapIdx])

        list(nRead=sum(nReadByTile),
             nFilter=sum(nFilterByTile),
             nClean=sum(nCleanByTile),
             nFilterClean=sum(nFilterCleanByTile),
             nMap=sum(nMapByTile),
             ##
             baseCall=rowSums(baseByCycle),
             qualityScore=rowSums(qualityByCycle),
             quality=quality,
             qualityFilter=qualityFilter,
             alignQuality=alignQuality,
             duplicates=list(
               tables=dupTables,
               mapTables=mapDupTables),
             perCycle=list(
               baseCall=baseByCycle,
               qualityScore=qualityByCycle),
             perTile=list(
               counts=data.frame(
                 Reads=nReadByTile,
                 Filtered=nFilterByTile,
                 Clean=nCleanByTile,
                 FilteredClean=nFilterCleanByTile,
                 Mapped=nMapByTile,
                 Tile=seq_along(nReadByTile),
                 Lane=pattern, row.names=NULL),
               qualityRead=tapply(qualityScore, df$tile, median),
               qualityFilter=tapply(qualityScore[filterIdx],
                 df$tile[filterIdx], median),
               qualityClean=tapply(
                 qualityScore[cleanIdx], df$tile[cleanIdx], median),
               qualityFilterClean=tapply(qualityScore[filterIdx & cleanIdx],
                 df$tile[filterIdx & cleanIdx], median),
               qualityMap=tapply(qualityScore[mapIdx], df$tile[mapIdx],
                 median)))
    }
    fls <- .file_names(dirPath, pattern)
    lst <- srapply(basename(fls), .lane, dirPath=dirPath, type=type)
    names(lst) <- basename(fls)

    .qualityDf <- function(lanes, part) {
        x <- lapply(lanes, function(elt) elt[[part]]$x)
        y <- lapply(lanes, function(elt) elt[[part]]$y)
        data.frame(Quality=unlist(x, use.names=FALSE),
                   Density=unlist(y, use.names=FALSE),
                   Lane=rep(names(lanes),
                     sapply(x, length, USE.NAMES=FALSE)))
    }

    alignQuality <- mapply(function(elt, lane) {
        data.frame(Quality=as.integer(names(elt[["alignQuality"]])),
                   Count=as.vector(elt[["alignQuality"]]),
                   Lane=lane)
    }, lst, names(lst), SIMPLIFY=FALSE, USE.NAMES=FALSE)

    .dupTop <- function(elt, lane) {
        data.frame(Sequence=names(elt$top),
                   Count=as.vector(elt$top),
                   Lane=lane, row.names=NULL)
    }
    .dupTable <- function(elt, lane) {
        data.frame(NUniqueSequences=elt$distribution$nUniqueSequences,
                   Count=elt$distribution$nReads,
                   Lane=lane)
    }
    dups <- lapply(lst, "[[", "duplicates")
    top <- mapply(.dupTop, lapply(dups, "[[", "tables"),
                  names(dups), SIMPLIFY=FALSE, USE.NAMES=FALSE)
    mapTop <- mapply(.dupTop, lapply(dups, "[[", "mapTables"),
                     names(dups), SIMPLIFY=FALSE, USE.NAMES=FALSE)
    tables <- mapply(.dupTable,
                     lapply(dups, "[[", "tables"),
                     names(dups), SIMPLIFY=FALSE, USE.NAMES=FALSE)
    mapTables <- mapply(.dupTable,
                        lapply(dups, "[[", "mapTables"),
                        names(dups), SIMPLIFY=FALSE, USE.NAMES=FALSE)
    duplicates <- list(top=do.call("rbind", top),
                       mapTop=do.call("rbind", mapTop),
                       tables=do.call("rbind", tables),
                       mapTables=do.call("rbind", mapTables))

    tileCounts <- lapply(lst, function(elt) elt$perTile[["counts"]])
    names(tileCounts) <- NULL

    list(counts=sapply(c("nRead", "nFilter", "nClean",
           "nFilterClean", "nMap"),
           function(elt) sapply(lst, "[[", elt)),
         baseCall=t(sapply(lst, "[[", "baseCall")),
         qualityRead=.qualityDf(lst, "quality"),
         qualityFilter=.qualityDf(lst, "qualityFilter"),
         alignQuality=do.call("rbind", alignQuality),
         duplicates=duplicates,
         perCycle=lapply(lst, "[[", "perCycle"),
         perTile=list(
           counts=do.call("rbind", tileCounts),
           quality=lapply(lst, function(elt) {
               as.vector(elt$perTile[names(elt$perTile) != "counts"])
           })))
}

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
