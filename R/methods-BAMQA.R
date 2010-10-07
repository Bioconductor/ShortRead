
.BAMQA <- function(x, ...)
{
    new("BAMQA", .srlist=x, ...)
}

.qa_BAM_lane <-
   function(dirPath, pattern, ...,  type="BAM", verbose=FALSE)
{
    if (verbose)
        message("qa 'BAM' pattern:", pattern)
    rpt <- readAligned(dirPath, pattern, type, ...)
	doc <- .qa_depthOfCoverage(rpt, pattern)
    # adapter contamination is computed in qa() 
    res <-  qa(rpt, pattern, ..., verbose=verbose)
    c(.srlist(res), 
		list(depthOfCoverage=doc)) 
}


.qa_BAM <-
    function(dirPath, pattern, type="BAM",..., 
			 param=list(ScanBamParam(
               	simpleCigar=TRUE,
               	reverseComplement=TRUE,
               	what=.readAligned_bamWhat())),
			 verbose=FALSE)
	
{
    fls <- .file_names(dirPath, pattern)
	if(!missing(param)) {
 		## make param into a list
        if (is.list(param) && !all(sapply(param, is, "ScanBamParam"))) {
            msg <-  paste("All elements of the param",
                          "argument must be ScanBamParm objects.",
                          collapse="")
            .throw(SRError("UserArgumentMismatch",msg))
        }
        else if (!is.list(param) && !is(param, "ScanBamParam")) {
            msg <-  paste("The param argument must be a",
                          "ScanBamParm object.",
                          collapse="")
            .throw(SRError("UserArgumentMismatch",msg))
		}
        else if (!is.list(param) && is(param, "ScanBamParam")) {
            param <- list(param)
		}

        if (length(param) > 1L && (length(param) != length(fls))) {
            msg <-  paste("The length of the param argument does not",
                          "match the number of BAM files to be read.",
                          " If more than one param is supplied there must",
                          "be one for each file.", collapse="")
            .throw(SRError("UserArgumentMismatch",msg))
        }

 		## FIXME: currently we only deal with cigars without indels
        if (any(sapply(param, bamSimpleCigar) != TRUE)) {
            msg <- paste("using 'TRUE' for 'bamSimpleCigar(param)'",
                "(skipping reads with I, D, H, S, or P in 'cigar')")
            .throw(SRWarn("UserArgumentMismatch", msg))
        }
        if (any(sapply(param, bamReverseComplement) != TRUE)) {
            msg <- "using 'TRUE' for 'bamReverseComplement(param)'"
            .throw(SRWarn("UserArgumentMismatch", msg))
        }
        if (any(sapply(param,
                FUN = function(x) {
                !setequal(bamWhat(x), .readAligned_bamWhat())
                }) == TRUE)) {
            msg <- sprintf("using '%s' for 'bamWhat(param)'",
                   paste(.readAligned_bamWhat(),
                   collapse="', '"))
            .throw(SRWarn("UserArgumentMismatch", msg))
        }

     	param <- sapply(param,
                    	FUN = function(param) {
                    	initialize(param, simpleCigar=TRUE,
                    	reverseComplement=TRUE,
                    	what=.readAligned_bamWhat())})

	}
        DF <- DataFrame(file=fls)
        DF$param <- param

	lst <- mapply(.qa_BAM_lane, pattern=basename(DF$file), 
				  param=DF$param, MoreArgs=list(dirPath=dirPath, ...), 
			      SIMPLIFY=FALSE, USE.NAMES=FALSE, 
				  verbose=verbose)

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
    qa <- x                             # mnemonic alias
    dir.create(dest, recursive=TRUE)
    fls <- c("0000-Header.html", "1000-Overview.html",
             "2000-RunSummary.html", "3000-ReadDistribution.html",
             "4000-CycleSpecific.html", "8000-DepthOfCoverage.html",
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
               .plotDepthOfCoverage(qa[["depthOfCoverage"]])),
			 ADAPTER_CONTAMINATION=hwrite(
			   .ppnCount(qa[["adapterContamination"]]),
			   border=NULL)
             )
    .report_html_do(dest, sections, values, ...)
}

setMethod(.report_html, "BAMQA", .report_html_BAMQA)

