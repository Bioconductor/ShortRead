.BAMQA <- function(x, ...)
{
    new("BAMQA", .srlist=x, ...)
}

.qa_BAM_lane <-
   function(dirPath, ..., verbose=FALSE)
{
    if (verbose)
        message("qa 'BAM' dirPath:", dirPath)
    rpt <- readAligned(dirPath, type="BAM", ...)
    fileName <- basename(dirPath)
    doc <- .qa_depthOfCoverage(rpt, fileName)
    res <-  .srlist(qa(rpt, fileName, ..., verbose=verbose))
    flag <- alignData(rpt)[["flag"]]
    res[["readCounts"]][,c("filter", "aligned")] <-
        c(sum(bamFlagTest(flag, "isValidVendorRead")),
          sum(!bamFlagTest(flag, "isUnmappedQuery")))
    c(res, list(depthOfCoverage=doc))
}

.qa_BAM <-
    function(dirPath, pattern, type="BAM", ...,
             param=ScanBamParam(simpleCigar=TRUE,
               reverseComplement=TRUE, what=.readAligned_bamWhat(FALSE)))
{
    fls <- .file_names(dirPath, pattern)
    lst <- srapply(fls, .qa_BAM_lane, ..., param=param, reduce=.reduce(1))
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
             adapterContamination=.bind(lst, "adapterContamination")
    )
    .BAMQA(lst)
}

.report_html_BAMQA <-
    function(x, dest, type, ...)
{
    qa <- .qa_sampleKey(x)
    dir.create(dest, recursive=TRUE)
    fls <- c("0000-Header.html", "1000-Overview.html",
             "2000-RunSummary.html", "3000-ReadDistribution.html",
             "4000-CycleSpecific.html", "8000-DepthOfCoverage.html",
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
             DEPTH_OF_COVERAGE_FIGURE=.html_img(
               dest, "depthOfCoverage",
               .plotDepthOfCoverage(qa[["depthOfCoverage"]])),
             ADAPTER_CONTAMINATION=hwrite(
             .df2a(qa[["adapterContamination"]]),
             border=0)
             )
    .report_html_do(dest, sections, values, ...)
}

setMethod(report_html, "BAMQA", .report_html_BAMQA)

