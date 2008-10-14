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

.qa_Solexa_tileStats <- function(dirPath, pattern, ...)
{
    .qa_Solexa_tileStats_tile <- function(dirPath, pattern, ...)
    {
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
    if (length(pattern)==0) pattern=".*_seq.txt$"
    lst <- srapply(list.files(dirPath, pattern, full.names=TRUE),
                   .qa_Solexa_tileStats_tile, dirPath=dirPath)
    .qa_lst_as_data_frame(lst)
}

.qa_SolexaExport_lane <-
    function(dirPath, pattern, ..., type="SolexaExport", verbose=FALSE)
{
    readLbls <- c("read", "aligned", "filtered")
    rpt <- readAligned(dirPath, pattern,..., type=type)
    df <- pData(alignData(rpt))
    filterIdx <- df$filtering=="Y"
    mapIdx <- !is.na(position(rpt))

    nbins <- max(df$tile)
    tiles <- seq_len(nbins)
    nReadByTile <- tabulate(df$tile, nbins)
    nFilterByTile <- tabulate(df$tile[filterIdx], nbins)
    nMapByTile <- tabulate(df$tile[mapIdx], nbins)

    qualityScore <- alphabetScore(quality(rpt)) / width(quality(rpt))
    qualityDf <- function(qscore)
    {
        d <- density(qscore)
        data.frame(quality=d$x,
                   density=d$y,
                   lane=pattern)
    }
    qualityScoreRead <- qualityDf(qualityScore)
    qualityScoreFiltered <- qualityDf(qualityScore[filterIdx])
    qualityScoreAligned <- qualityDf(qualityScore[mapIdx])

    abc <- alphabetByCycle(rpt)
    baseQuality <- apply(abc, 2, sum)

    alignQuality <- table(quality(alignQuality(rpt))[mapIdx])

    tablesRead <- tables(sread(rpt))
    tablesFiltered <- tables(sread(rpt)[filterIdx])
    tablesAligned <- tables(sread(rpt)[mapIdx])
    frequentSequences <-
        data.frame(sequence=c(
                     names(tablesRead$top),
                     names(tablesFiltered$top),
                     names(tablesAligned$top)),
                   count=c(
                     as.integer(tablesRead$top),
                     as.integer(tablesFiltered$top),
                     as.integer(tablesAligned$top)),
                   type=rep(
                     readLbls,
                     c(length(tablesRead$top),
                       length(tablesFiltered$top),
                       length(tablesAligned$top))),
                   lane=pattern)
    sequenceDistribution <-
        cbind(rbind(tablesRead$distribution,
                    tablesFiltered$distribution,
                    tablesAligned$distribution),
              type=rep(
                readLbls,
                c(nrow(tablesRead$distribution),
                  nrow(tablesFiltered$distribution),
                  nrow(tablesAligned$distribution))),
              lane=pattern)

    perCycleBaseCall <- local({
        abc <- apply(abc, c(1, 3), sum)
        df <- data.frame(Cycle=col(abc, as.factor=TRUE),
                         Base=row(abc, as.factor=TRUE),
                         Count=as.vector(abc),
                         lane=pattern)
        df[df$Count != 0,]
    })
    perCycleQuality <- local({
        abc <- apply(abc, 2:3, sum)
        q <- row(abc, as.factor=TRUE)
        q0 <- 1 + 32 * is(quality(rpt), "SFastqQuality")
        df <- data.frame(Cycle=col(abc, as.factor=TRUE),
                         Quality=q,
                         Score=as.numeric(q)-q0,
                         Count=as.vector(abc),
                         lane=pattern)
        df[df$Count != 0, ]
    })

    list(readCounts=data.frame(
           read=sum(nReadByTile),
           filtered=sum(nFilterByTile),
           aligned=sum(nMapByTile),
           row.names=pattern),
         baseCalls=local({
             n <- apply(abc, 1, sum)
             data.frame(A=n["A"], C=n["C"], G=n["G"], T=n["T"],
                        N=n["N"], row.names=pattern)
         }),
         readQualityScore=cbind(
           rbind(qualityScoreRead,
                 qualityScoreFiltered,
                 qualityScoreAligned),
           type=rep(
             readLbls,
             c(nrow(qualityScoreRead),
               nrow(qualityScoreFiltered),
               nrow(qualityScoreAligned)))),
         baseQuality=data.frame(
           score=as.vector(names(baseQuality)),
           count=as.vector(baseQuality),
           lane=pattern, row.names=NULL),
         alignQuality=data.frame(
           score=as.numeric(names(alignQuality)),
           count=as.vector(alignQuality),
           lane=pattern, row.names=NULL),
         frequentSequences=frequentSequences,
         sequenceDistribution=sequenceDistribution,

         perCycle=list(
           baseCall=perCycleBaseCall,
           quality=perCycleQuality),

         perTile=list(
           readCounts=data.frame(
             count=c(nReadByTile, nFilterByTile, nMapByTile),
             type=rep(
               readLbls,
               c(length(nReadByTile), length(nFilterByTile),
                 length(nMapByTile))),
             tile=rep(tiles, 3),
             lane=pattern, row.names=NULL),
           medianReadQualityScore=local({
               tidx <- as.character(tiles)
               data.frame(score=c(
                            tapply(qualityScore,
                                   df$tile, median)[tidx],
                            tapply(qualityScore[filterIdx],
                                   df$tile[filterIdx], median)[tidx],
                            tapply(qualityScore[mapIdx],
                                   df$tile[mapIdx], median)[tidx]),
                          type=rep(readLbls, each=length(tidx)),
                          tile=as.integer(tidx),
                          lane=pattern, row.names=NULL)
           }))
         )
}

.qa_SolexaExport <- function(dirPath, pattern, type="SolexaExport", ...,
                              verbose=FALSE) {
    fls <- .file_names(dirPath, pattern)
    lst <- srapply(basename(fls), .qa_SolexaExport_lane,
                   dirPath=dirPath, type=type,
                   verbose=verbose)
    names(lst) <- basename(fls)

    ## collapse into data frames
    bind <- function(lst, elt)
        do.call("rbind",
                subListExtract(lst, elt, keep.names=FALSE))
    lst <-
        list(readCounts=bind(lst, "readCounts"),
             baseCalls=bind(lst, "baseCalls"),
             readQualityScore=bind(lst, "readQualityScore"),
             baseQuality=bind(lst, "baseQuality"),
             alignQuality=bind(lst, "alignQuality"),
             frequentSequences=bind(lst, "frequentSequences"),
             sequenceDistribution=bind(lst, "sequenceDistribution"),
             perCycle=local({
                 lst <- subListExtract(lst, "perCycle")
                 list(baseCall=bind(lst, "baseCall"),
                      quality=bind(lst, "quality"))
             }),
             perTile=local({
                 lst <- subListExtract(lst, "perTile")
                 list(readCounts=bind(lst, "readCounts"),
                      medianReadQualityScore=bind(
                        lst, "medianReadQualityScore"))
             }))
    .SolexaExportQA(lst)
}

## .qa_SolexaSet <- function(dirPath, pattern=character(0), ...)
## {
##     dirPath <- analysisPath(dirPath)
##     if (missing(pattern))
##         pattern <- ".*_export.txt$"
##     callGeneric(dirPath, pattern, type="SolexaExport", ...)
## }

## setMethod("qa", "SolexaSet", .qa_solexa_export)

## .report_SolexaSet <- function(x, run=1, ..., qaFile=tempfile(),
##                               dest=tempfile(), type="pdf" )
## {
##     report(qa(x, run=run))
## }

## setMethod("report", "SolexaSet", .report_SolexaSet)


## alignment

.readAligned_SolexaSet <- function(dirPath,
                                   pattern=".*_export.txt$",
                                   run, ...)
{
    dirPath <- analysisPath(solexaPath(dirPath))[run]
    .readAligned_character(dirPath, pattern, type="SolexaExport", ...)
}

setMethod(readAligned, "SolexaSet", .readAligned_SolexaSet)

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
