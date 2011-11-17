.SolexaExportQA <- function(x, ...) 
{
    new("SolexaExportQA", .srlist=x, ...)
}

## qa

.qa_lst_as_data_frame <- function(lst) {
    if (length(lst)==0) return(data.frame())
    nms <- names(lst[[1]])
    sublst <- sapply(nms, function(nm) {
        subListExtract(lst, nm, simplify=TRUE)
    })
    names(sublst) <- nms
    do.call(data.frame, sublst)
}

.qa_Solexa_tileStats <- function(dirPath, pattern, ...)
{
    .qa_Solexa_tileStats_tile <- function(dirPath, pattern, ...)
    {
        lane <- as.numeric(sub("s_([0-9]+)_.*", "\\1", pattern))
        tile <- as.numeric(sub("s_[0-9]+_([0-9]+)_.*", "\\1",
                               pattern))
        dna <- readXStringColumns(dirPath, pattern,
                                  colClasses=list(NULL, NULL, NULL,
                                    NULL, "DNAString"))[[1]]
        list(lane=lane, tile=tile,
             slane=(lane-1)*3+trunc((tile-1)/100)+1,
             stile=1+pmin((tile-1)%%200, (200-tile)%%200),
             nReads=length(dna),
             nClean=sum(alphabetFrequency(dna, baseOnly=TRUE)[,"other"]==0))
    }
    if (length(pattern)==0) pattern=".*_seq.txt$"
    lst <- srapply(list.files(dirPath, pattern, full.names=TRUE),
                   .qa_Solexa_tileStats_tile, dirPath=dirPath, ...)
    .qa_lst_as_data_frame(lst)
}

.qa_SolexaExport_lane <-
    function(dirPath, pattern, ..., type="SolexaExport",
             verbose=FALSE)
{
    if (verbose)
        message("qa 'SolexaExport' pattern:", pattern)
    readLbls <- c("read", "filtered", "aligned")
    rpt <- readAligned(dirPath, pattern, type, ...)
    doc <- .qa_depthOfCoverage(rpt, pattern)
    ac <- .qa_adapterContamination(rpt, pattern, ...)
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
    perCycleBaseCall <- .qa_perCycleBaseCall(abc, pattern)
    perCycleQuality <- .qa_perCycleQuality(abc, quality(rpt), pattern)
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
           })),
         depthOfCoverage=doc,
         adapterContamination=ac)
}

.qa_SolexaExport <-
    function(dirPath, pattern, type="SolexaExport", ...,
             verbose=FALSE)
{
    fls <- .file_names(dirPath, pattern)
    lst <- srapply(basename(fls), .qa_SolexaExport_lane,
                   dirPath=dirPath, type=type, ...,
                   reduce=.reduce(1), verbose=verbose, 
                   USE.NAMES=TRUE)

    ## collapse into data frames
    lst <-
        list(readCounts=.bind(lst, "readCounts"),
             baseCalls=.bind(lst, "baseCalls"),
             readQualityScore=.bind(lst, "readQualityScore"),
             baseQuality=.bind(lst, "baseQuality"),
             alignQuality=.bind(lst, "alignQuality"),
             frequentSequences=.bind(lst, "frequentSequences"),
             sequenceDistribution=.bind(lst, "sequenceDistribution"),
             perCycle=local({
                 lst <- subListExtract(lst, "perCycle")
                 list(baseCall=.bind(lst, "baseCall"),
                      quality=.bind(lst, "quality"))
             }),
             perTile=local({
                 lst <- subListExtract(lst, "perTile")
                 list(readCounts=.bind(lst, "readCounts"),
                      medianReadQualityScore=.bind(
                        lst, "medianReadQualityScore"))
             }),
             depthOfCoverage=.bind(lst, "depthOfCoverage"),
             adapterContamination=.bind(lst, "adapterContamination"))
    .SolexaExportQA(lst)
}

## report

setMethod(.report_pdf, "SolexaExportQA",
          function(x, dest, type, ...)
{
    qa <- x                             # mnemonic alias
    to <- tempfile()
    save(qa, file=to)
    res <- callGeneric(to, dest, type, ...)
    unlink(to)
    res
})

setMethod(report_html, "SolexaExportQA",
          function (x, dest, type, ...)
{
    qa <- .qa_sampleKey(x)
    dir.create(dest, recursive=TRUE)
    fls <- c("0000-Header.html", "1000-Overview.html",
             "2000-RunSummary.html", "3000-ReadDistribution.html",
             "4000-CycleSpecific.html", "5000-PerTile.html",
             "6000-Alignment.html", "8000-DepthOfCoverage.html",
             "9000-AdapterContamination.html", "9999-Footer.html")
    sections <- system.file("template", fls, package="ShortRead")
    perCycle <- qa[["perCycle"]]
    perTile <- qa[["perTile"]]
    readCnt <- perTile[["readCounts"]]
    values <-
        list(SAMPLE_KEY=hwrite(qa[["keyValue"]], border=0),
             PPN_COUNT=.html_img(
               dest, "readCount", .plotReadCount(qa)),
             PPN_COUNT_TBL=hwrite(
               .ppnCount(qa[["readCounts"]]),
               border=0),
             BASE_CALL_COUNT=.html_img(
               dest, "baseCalls", .plotNucleotideCount(qa)),
             READ_QUALITY_FIGURE=.htmlReadQuality(
               dest, "readQuality", qa),
             READ_OCCURRENCES_FIGURE=.htmlReadOccur(
               dest, "readOccurences", qa),
             FREQUENT_SEQUENCES_READ=hwrite(
               .freqSequences(qa, "read"),
               border=0),
             FREQUENT_SEQUENCES_FILTERED=hwrite(
               .freqSequences(qa, "filtered"),
               border=0),
             FREQUENT_SEQUENCES_ALIGNED=hwrite(
               .freqSequences(qa, "aligned"),
               border=0),
             CYCLE_BASE_CALL_FIGURE=.html_img(
               dest, "perCycleBaseCall",
               .plotCycleBaseCall(perCycle$baseCall)),
             CYCLE_QUALITY_FIGURE=.html_img(
               dest, "perCycleQuality",
               .plotCycleQuality(perCycle$quality)),
             PER_TILE_HISTOGRAM=local({
                 cnts <- readCnt[readCnt$type=="read", "count"]
                 hist <- histogram(cnts, breaks=40,
                                   xlab="Reads per tile",
                                   panel=function(x, ...) {
                                       panel.abline(v=quantile(x, .1),
                                                    col="red", lty=2)
                                       panel.histogram(x, ...)
                                   }, col="white")
                 .html_img(dest, "perTileHistogram", hist)
             }),
             PER_TILE_COUNT_FIGURE=.html_img(
               dest, "perTileCount",
               .plotTileCounts(readCnt[readCnt$type=="read",])),
             PER_TILE_QUALITY_FIGURE=local({
                 qscore <- perTile[["medianReadQualityScore"]]
                 score <- qscore[qscore$type=="read",]
                 .html_img(dest, "perTileQuality",
                           .plotTileQualityScore(score))
             }),
             ALIGN_QUALITY_FIGURE=.html_img(
               dest, "alignmentQuality",
               .plotAlignQuality(qa[["alignQuality"]])),
             DEPTH_OF_COVERAGE_FIGURE=.html_img(
               dest, "depthOfCoverage",
               .plotDepthOfCoverage(qa[["depthOfCoverage"]])),
             ADAPTER_CONTAMINATION=hwrite(
               .df2a(qa[["adapterContamination"]]),
               border=0))
    .report_html_do(dest, sections, values, ...)
})
