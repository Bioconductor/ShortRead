.MAQMapQA <- function(x, ...)
{
    new("MAQMapQA", .srlist=x, ...)
}

.maq_reverse <- function(aln)
{
    plus <- strand(aln) == "+"
    new("AlignedRead",
        append(aln[plus], aln[!plus]),
        sread=append(
          sread(aln)[plus],
          reverseComplement(sread(aln)[!plus])),
        quality=append(
          quality(aln)[plus],
          FastqQuality(reverse(quality(quality(aln)[!plus])))))
}

.qa_MAQMap_lane <-
    function(dirPath, pattern, type, ..., verbose=FALSE) 
{
    if (verbose)
        message("qa '", type, "' pattern: ", pattern, sep="")
    rpt <- .maq_reverse(readAligned(dirPath, pattern, type, ...))
    alf <- alphabetFrequency(sread(rpt), baseOnly=TRUE,collapse=TRUE)
    bqtbl <- alphabetFrequency(quality(rpt), collapse=TRUE)
    rqs <- .qa_qdensity(quality(rpt))
    freqtbl <- tables(sread(rpt))
    abc <- alphabetByCycle(rpt)
    perCycleBaseCall <- .qa_perCycleBaseCall(abc, pattern)
    perCycleQuality <- .qa_perCycleQuality(abc, quality(rpt), pattern)
    aqtbl <- table(quality(alignQuality(rpt)))
    list(readCounts=data.frame(
           read=NA, filter=NA, aligned=length(rpt),
           row.names=pattern),
         baseCalls=data.frame(
           A=alf[["A"]], C=alf[["C"]], G=alf[["G"]], T=alf[["T"]],
           N=alf[["other"]], row.names=pattern),
         readQualityScore=data.frame(
           quality=rqs$x,
           density=rqs$y,
           lane=pattern,
           type="aligned", row.names=NULL),
         baseQuality=data.frame(
           score=names(bqtbl),
           count=as.vector(bqtbl),
           lane=pattern, row.names=NULL),
         alignQuality=data.frame(
           score=as.numeric(names(aqtbl)),
           count=as.vector(aqtbl),
           lane=pattern, row.names=NULL),
         frequentSequences=data.frame(
           sequence=names(freqtbl$top),
           count=as.integer(freqtbl$top),
           type="aligned",
           lane=pattern, row.names=NULL),
         sequenceDistribution=cbind(
           freqtbl$distribution,
           type="aligned",
           lane=pattern, row.names=NULL),
         perCycle=list(
           baseCall=perCycleBaseCall,
           quality=perCycleQuality),
         perTile=list(
           readCounts=data.frame(
             count=integer(0), type=character(0),
             tile=integer(0), lane=character(0)),
           medianReadQualityScore=data.frame(
             score=integer(), type=character(), tile=integer(),
             lane=integer()))
         )
}

.qa_MAQMap <-
    function(dirPath, pattern, type, ..., verbose=FALSE)
{
    fls <- .file_names(dirPath, pattern)
    lst <- srapply(basename(fls), .qa_MAQMap_lane,
                   dirPath=dirPath, type=type, ...,
                   reduce=.reduce(1), verbose=verbose, USE.NAMES=TRUE)
    
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
             }))
    .MAQMapQA(lst)
}

setMethod(report_html, "MAQMapQA",
          function(x, dest, type, ...)
{
    qa <- .qa_sampleKey(x)
    dir.create(dest, recursive=TRUE)
    fls <- c("0000-Header.html", "1000-Overview.html",
             "2000-RunSummary.html", "3000-ReadDistribution.html",
             "4000-CycleSpecific.html",
             "6000-Alignment.html",
             "9999-Footer.html")
    sections <- system.file("template", fls, package="ShortRead")
    perCycle <- qa[["perCycle"]]
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
               dest, "readQuality", qa, "aligned"),
             READ_OCCURRENCES_FIGURE=.htmlReadOccur(
               dest, "readOccurences", qa, "aligned"),
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
             ALIGN_QUALITY_FIGURE=.html_img(
               dest, "alignmentQuality",
               .plotAlignQuality(qa[["alignQuality"]]))
             )
    .report_html_do(dest, sections, values, ...)
})
          
