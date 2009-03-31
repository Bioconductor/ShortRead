.FastqQA <- function(x, ...)
{
    new("FastqQA", .srlist=x, ...)
}

.qa_fastq_lane <-
    function(dirPath, pattern, ..., type="fastq", verbose=FALSE)
{
    if (verbose)
        message("qa 'fastq' pattern:", pattern)
    rpt <- readFastq(dirPath, pattern)

    alf <- alphabetFrequency(sread(rpt), baseOnly=TRUE, collapse=TRUE)
    bqtbl <- alphabetFrequency(quality(rpt), collapse=TRUE)
    rqs <- local({
        qscore <- alphabetScore(quality(rpt)) / width(quality(rpt))
        density(qscore)
    })
    freqtbl <- tables(sread(rpt))
    abc <- alphabetByCycle(rpt)
    perCycleBaseCall <- local({
        abc <- apply(abc, c(1, 3), sum)
        df <- data.frame(Cycle=as.integer(colnames(abc)[col(abc)]),
                         Base=factor(rownames(abc)[row(abc)]),
                         Count=as.vector(abc),
                         lane=pattern)
        df[df$Count != 0,]
    })
    perCycleQuality <- local({
        abc <- apply(abc, 2:3, sum)
        q <- factor(rownames(abc)[row(abc)])
        q0 <- 1 + 32 * is(quality(rpt), "SFastqQuality")
        df <- data.frame(Cycle=as.integer(colnames(abc)[col(abc)]),
                         Quality=q,
                         Score=as.numeric(q)-q0,
                         Count=as.vector(abc),
                         lane=pattern)
        df[df$Count != 0, ]
    })
    list(readCounts=data.frame(
           read=length(rpt), filter=NA, aligned=NA,
           row.names=pattern),
         baseCalls=data.frame(
           A=alf[["A"]], C=alf[["C"]], G=alf[["G"]], T=alf[["T"]],
           N=alf[["other"]], row.names=pattern),
         readQualityScore=data.frame(
           quality=rqs$x,
           density=rqs$y,
           lane=pattern,
           type="read"),
         baseQuality=data.frame(
           score=names(bqtbl),
           count=as.vector(bqtbl),
           lane=pattern),
         alignQuality=data.frame(
           score=as.numeric(NA),
           count=as.numeric(NA),
           lane=pattern, row.names=NULL),
         frequentSequences=data.frame(
           sequence=names(freqtbl$top),
           count=as.integer(freqtbl$top),
           type="read",
           lane=pattern),
         sequenceDistribution=cbind(
           freqtbl$distribution,
           type="read",
           lane=pattern),
         perCycle=list(
           baseCall=perCycleBaseCall,
           quality=perCycleQuality),
         perTile=list(
           readCounts=data.frame(
             count=integer(0), type=character(0),
             tile=integer(0), lane=character(0)),
           medianReadQualityScore=data.frame(
             score=integer(), type=character(), tile=integer(),
             lane=integer(), row.names=NULL))
         )
}

.qa_fastq <-
    function(dirPath, pattern, type="fastq", ...,
             verbose=FALSE) 
{
    fls <- .file_names(dirPath, pattern)
    lst <- srapply(basename(fls), .qa_fastq_lane,
                   dirPath=dirPath, type=type,
                   verbose=verbose)
    names(lst) <- basename(fls)
    bind <- function(lst, elt)
        do.call(rbind,
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
    .FastqQA(lst)
}

setMethod(.report_html, "FastqQA",
    function (x, dest, type, ...)
{
    qa <- x                             # mnemonic alias
    dir.create(dest, recursive=TRUE)
    fls <- c("0000-Header.html", "1000-Overview.html",
             "2000-RunSummary.html", "3000-ReadDistribution.html",
             "4000-CycleSpecific.html", "9999-Footer.html")
    sections <- system.file("template", fls, package="ShortRead")
    perCycle <- qa[["perCycle"]]
    values <-
        list(PPN_COUNT=hwrite(
               .ppnCount(qa[["readCounts"]]),
               border=NULL),
             BASE_CALL_COUNT=hwrite(
               .df2a(qa[["baseCalls"]] / rowSums(qa[["baseCalls"]])),
               border=NULL),
             READ_QUALITY_FIGURE=.htmlReadQuality(
               dest, "readQuality", qa),
             READ_OCCURRENCES_FIGURE=.htmlReadOccur(
               dest, "readOccurences", qa),
             FREQUENT_SEQUENCES_READ=hwrite(
               .freqSequences(qa, "read"),
               border=NULL),
             FREQUENT_SEQUENCES_FILTERED=.html_NA(),
             FREQUENT_SEQUENCES_ALIGNED=.html_NA(),
             CYCLE_BASE_CALL_FIGURE=.html_img(
               dest, "perCycleBaseCall",
               .plotCycleBaseCall(perCycle$baseCall)),
             CYCLE_QUALITY_FIGURE=.html_img(
               dest, "perCycleQuality",
               .plotCycleQuality(perCycle$quality))
             )
    .report_html_do(dest, sections, values, ...)
})
