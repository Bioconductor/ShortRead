.SolexaRealignQA <- function(x, ...) 
{
    new("SolexaRealignQA", .srlist=x, ...)
}

## qa

.qa_SolexaRealign_lane <-
    function(dirPath, pattern, ..., type="SolexaRealign",
             verbose=FALSE)
{
    if (verbose)
        message("qa 'SolexaRealign' pattern:", pattern)
    readLbls <- c("read", "aligned")
    aln <- readAligned(dirPath, pattern, type=type, ...)
	doc <- .qa_depthOfCoverage(aln, pattern)
    ac <- .qa_adapterContamination(aln, pattern, ...)
    df <- pData(alignData(aln))

    mapIdx <- alignData(aln)[["nMatch"]] == 1L

    alf <- alphabetFrequency(sread(aln), baseOnly=TRUE, collapse=TRUE)
    abc <- alphabetByCycle(aln)

    alignQuality <- table(quality(alignQuality(aln))[mapIdx])

    tablesRead <- tables(sread(aln))
    tablesAligned <- tables(sread(aln)[mapIdx])
    frequentSequences <-
        data.frame(sequence=c(
                     names(tablesRead$top),
                     names(tablesAligned$top)),
                   count=c(
                     as.integer(tablesRead$top),
                     as.integer(tablesAligned$top)),
                   type=rep(
                     readLbls,
                     c(length(tablesRead$top),
                       length(tablesAligned$top))),
                   lane=pattern)
    sequenceDistribution <-
        cbind(rbind(tablesRead$distribution,
                    tablesAligned$distribution),
              type=rep(
                readLbls,
                c(nrow(tablesRead$distribution),
                  nrow(tablesAligned$distribution))),
              lane=pattern)

    perCycleBaseCall <- .qa_perCycleBaseCall(abc, pattern)
    perCycleQuality <- .qa_perCycleQuality()
    malntbl <- table(alignData(aln)[["nMatch"]])

    list(readCounts=data.frame(
           read=length(aln), filtered=NA, aligned=sum(mapIdx),
           row.names=pattern),
         baseCalls=data.frame(
           A=alf[["A"]], C=alf[["C"]], G=alf[["G"]], T=alf[["T"]],
           N=alf[["other"]], row.names=pattern),
         readQualityScore=data.frame(
           score=numeric(0),
           type=factor(character(0), levels=readLbls)),
         baseQuality=data.frame(
           score=numeric(0), count=integer(0), lane=character(0),
           row.names=NULL),
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
             count=integer(0), type=character(0),
             tile=integer(0), lane=character(0)),
           medianReadQualityScore=data.frame(
             score=integer(), type=character(), tile=integer(),
             lane=integer(), row.names=NULL)),

         multipleAlignment=data.frame(
           Count=as.vector(malntbl),
           Matches=as.integer(names(malntbl)),
           lane=pattern, row.names=NULL),

		 depthOfCoverage=doc,
		 adapterContamination=ac
         )
}

.qa_SolexaRealign <-
    function(dirPath, pattern, type="SolexaRealign", ...,
			 verbose=FALSE)
{
    fls <- .file_names(dirPath, pattern)
    lst <- srapply(basename(fls), .qa_SolexaRealign_lane,
                   dirPath=dirPath, type=type, ...,
				   reduce=.reduce(1), verbose=verbose, 
				   USE.NAMES=TRUE)

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
             multipleAlignment=.bind(lst, "multipleAlignment"),
             depthOfCoverage=.bind(lst, "depthOfCoverage"),
             adapterContamination=.bind(lst, "adapterContamination"))
    .SolexaRealignQA(lst)
}

## report

setMethod(report_html, "SolexaRealignQA",
          function (x, dest, type, ...)
{
    qa <- x                             # mnemonic alias
    dir.create(dest, recursive=TRUE)
    fls <- c("0000-Header.html", "1000-Overview.html",
             "1100-Overview-SolexaRealign.html",
             "2000-RunSummary.html", "3000-ReadDistribution.html",
             "4000-CycleSpecific.html", "6000-Alignment.html",
             "7000-MultipleAlignment.html", "8000-DepthOfCoverage.html",
             "9000-AdapterContamination.html", "9999-Footer.html")
    sections <- system.file("template", fls, package="ShortRead")
    perCycle <- qa[["perCycle"]]
    values <-
        list(PPN_COUNT=hwrite(
               .ppnCount(qa[["readCounts"]]),
               border=NULL),
             BASE_CALL_COUNT=hwrite(
               .df2a(qa[["baseCalls"]] / rowSums(qa[["baseCalls"]])),
               border=NULL),
             READ_QUALITY_FIGURE=.html_NA(),
             READ_OCCURRENCES_FIGURE=.htmlReadOccur(
               dest, "readOccurences", qa),
             FREQUENT_SEQUENCES_READ=hwrite(
               .freqSequences(qa, "read"),
               border=NULL),
             FREQUENT_SEQUENCES_FILTERED=.html_NA(),
             FREQUENT_SEQUENCES_ALIGNED=hwrite(
               .freqSequences(qa, "aligned"),
               border=NULL),
             CYCLE_BASE_CALL_FIGURE=.html_img(
               dest, "perCycleBaseCall",
               .plotCycleBaseCall(perCycle$baseCall)),
             CYCLE_QUALITY_FIGURE=.html_NA(),
             ALIGN_QUALITY_FIGURE=.html_img(
               dest, "alignmentQuality",
               .plotAlignQuality(qa[["alignQuality"]])),
             MULTIPLE_ALIGNMENT_COUNT_FIGURE=.html_img(
               dest, "multipleAlignmentCount",
               .plotMultipleAlignmentCount(qa[["multipleAlignment"]])),
             DEPTH_OF_COVERAGE_FIGURE=.html_img(
               dest, "depthOfCoverage",
               .plotDepthOfCoverage(qa[["depthOfCoverage"]])),
             ADAPTER_CONTAMINATION=hwrite(
               .ppnCount(qa[["adapterContamination"]]),
               border=NULL))

    .report_html_do(dest, sections, values, ...)
})
