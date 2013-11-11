.ShortReadQQA <- function(x, ...)
{
    new("ShortReadQQA", .srlist=x, ...)
}

.FastqQA <- function(x, ...)
{
    new("FastqQA", .srlist=x, ...)
}

.qa_ShortReadQ <-
    function(dirPath, lane, ..., verbose=FALSE)
{
    if (missing(lane))
        .throw(SRError("UserArgumentMismatch",
                       "'%s' must be '%s'", "lane", "character(1)"))
    obj <- dirPath 
    alf <- .qa_alphabetFrequency(sread(obj), baseOnly=TRUE, collapse=TRUE)
    bqtbl <- .qa_alphabetFrequency(quality(obj), collapse=TRUE)
    rqs <- .qa_qdensity(quality(obj))
    freqtbl <- tables(sread(obj))
    abc <- alphabetByCycle(obj)
    ac <- .qa_adapterContamination(obj, lane, ...)
    perCycleBaseCall <- .qa_perCycleBaseCall(abc, lane)
    perCycleQuality <- .qa_perCycleQuality(abc, quality(obj), lane)
    lst <- list(readCounts=data.frame(
           read=length(obj), filter=NA, aligned=NA,
           row.names=lane),
         baseCalls=data.frame(
           A=alf[["A"]], C=alf[["C"]], G=alf[["G"]], T=alf[["T"]],
           N=alf[["other"]], row.names=lane),
         readQualityScore=data.frame(
           quality=rqs$x,
           density=rqs$y,
           lane=lane,
           type="read"),
         baseQuality=data.frame(
           score=names(bqtbl),
           count=as.vector(bqtbl),
           lane=lane),
         alignQuality=data.frame(
           score=as.numeric(NA),
           count=as.numeric(NA),
           lane=lane, row.names=NULL),
         frequentSequences=data.frame(
           sequence=names(freqtbl$top),
           count=as.integer(freqtbl$top),
           type="read",
           lane=lane),
         sequenceDistribution=cbind(
           freqtbl$distribution,
           type="read",
           lane=lane),
         perCycle=list(
           baseCall=perCycleBaseCall,
           quality=perCycleQuality),
         perTile=list(
           readCounts=data.frame(
             count=integer(0), type=character(0),
             tile=integer(0), lane=character(0)),
           medianReadQualityScore=data.frame(
             score=integer(), type=character(), tile=integer(),
             lane=integer(), row.names=NULL)),
         adapterContamination=ac

         )

    .ShortReadQQA(lst)
}

setMethod(qa, "ShortReadQ", .qa_ShortReadQ)

.qa_fastq_lane <-
    function(dirPath, ..., sample=TRUE, type="fastq", 
        verbose=FALSE)
{
    if (verbose)
        message("qa 'fastq' pattern:", pattern)
    if (sample) {
        samp <- FastqSampler(dirPath, ...)
        qa <- qa(yield(samp), basename(dirPath), ..., verbose=verbose)
        close(samp)
        elts <- .srlist(qa)
        elts$readCounts$read <- samp$status()[["total"]]
        initialize(qa, .srlist=elts)
    } else {
        fq <-readFastq(dirPath, pattern, ...)
        qa(fq, pattern, ..., verbose=verbose)
    }
}

.qa_fastq <-
    function(dirPath, pattern, type="fastq", ...,
        verbose=FALSE) 
{
    fls <- .file_names(dirPath, pattern)
    lst <- bplapply(fls, .qa_fastq_lane, type=type, ...,
                    verbose=verbose)
    lst <- do.call(rbind, lst)
    .FastqQA(.srlist(lst))              # re-cast
}

.report_html_ShortReadQA <-             # or FastqQA
    function(x, dest, type, ...)
{
    qa <- .qa_sampleKey(x)
    dir.create(dest, recursive=TRUE)
    fls <- c("0000-Header.html", "1000-Overview.html",
             "2000-RunSummary.html", "3000-ReadDistribution.html",
             "4000-CycleSpecific.html", 
             "9000-AdapterContamination.html", "9999-Footer.html")
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
               dest, "readQuality", qa),
             READ_OCCURRENCES_FIGURE=.htmlReadOccur(
               dest, "readOccurences", qa),
             FREQUENT_SEQUENCES_READ=hwrite(
               .freqSequences(qa, "read"),
               border=0),
             FREQUENT_SEQUENCES_FILTERED=.html_NA(),
             FREQUENT_SEQUENCES_ALIGNED=.html_NA(),
             CYCLE_BASE_CALL_FIGURE=.html_img(
               dest, "perCycleBaseCall",
               .plotCycleBaseCall(perCycle$baseCall)),
             CYCLE_QUALITY_FIGURE=.html_img(
               dest, "perCycleQuality",
               .plotCycleQuality(perCycle$quality)),
             ADAPTER_CONTAMINATION=hwrite(
               .df2a(qa[["adapterContamination"]]),
               border=0)

             )
    .report_html_do(dest, sections, values, ...)
}

setMethod(report_html, "ShortReadQQA", .report_html_ShortReadQA)
setMethod(report_html, "FastqQA", .report_html_ShortReadQA)

