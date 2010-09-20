
.BAMQA <- function(x, ...)
{
    new("BAMQA", .srlist=x, ...)
}



.qa_BAM_lane <-
   function(dirPath, pattern, ..., type="Bowtie", verbose=FALSE)
{
    if (verbose)
        message("qa 'BAM' pattern:", pattern)
    rpt <- readAligned(dirPath, pattern, type, ...)
    res <-  qa(rpt, pattern, verbose=verbose)
	doc <- .qa_depthOfCoverage(rpt[occurrenceFilter(withSread=FALSE)(rpt)],
                               pattern)

    list(readCounts=res[["readCounts"]],
         baseCalls=res[["baseCalls"]],
         readQualityScore=res[["readQualityScore"]],
         baseQuality=res[["baseQuality"]],
         alignQuality=res[["alignQuality"]],
         frequentSequences=res[["frequentSequences"]],
         sequenceDistribution=res[["sequenceDistribution"]],
         perCycle=res[["perCycle"]],
         perTile=res[["perTile"]],
		 depthOfCoverage=doc
         )
}


.qa_BAM <-
    function(dirPath, pattern, type="BAM",
    	..., verbose=FALSE)
{
    fls <- .file_names(dirPath, pattern)
 	lst <- srapply(basename(fls), .qa_BAM_lane,
                   dirPath=dirPath, type=type, ...,
                   reduce=.reduce(1), verbose=verbose, USE.NAMES=TRUE)
    #lst <- do.call(rbind, lst)
    #.BAMQA(.srlist(lst))              # re-cast
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
             }),
             depthOfCoverage=bind(lst, "depthOfCoverage")
		)
    .BAMQA(lst)
}



.report_html_BAMQA <-          
    function(x, dest, type, ...)
{
    qa <- x                             # mnemonic alias
    dir.create(dest, recursive=TRUE)
    fls <- c("0000-Header.html", "1000-Overview.html",
             "2000-RunSummary.html", "3000-ReadDistribution.html",
             "4000-CycleSpecific.html", "8000-DepthOfCoverage.html",
			 "9999-Footer.html")
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
               .plotCycleQuality(perCycle$quality)),
             DEPTH_OF_COVERAGE_FIGURE=.html_img(
               dest, "depthOfCoverage",
               .plotDepthOfCoverage(qa[["depthOfCoverage"]]))
             )
    .report_html_do(dest, sections, values, ...)
}

setMethod(.report_html, "BAMQA", .report_html_BAMQA)

