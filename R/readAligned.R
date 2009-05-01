.read_csv_portion <- function(dirPath, pattern, colClasses, ...) {
    ## visit several files, then collapse
    files <- .file_names(dirPath, pattern)
    lsts <- lapply(files, function(fl, ...) {
        tryCatch({
            read.csv(fl, ...)
        }, error=function(err) {
            read.csv(gzfile(fl), ...)
        })
    }, ..., colClasses=colClasses, stringsAsFactors=FALSE)
    cclasses <- colClasses[!sapply(colClasses, is.null)]
    lst <- lapply(seq_along(names(cclasses)),
                  function(idx) unlist(lapply(lsts, "[[", idx)))
    names(lst) <- names(cclasses)
    lst
}

.readAligned_SolexaAlign <-
    function(dirPath, pattern=character(0), ...,
             quote="", sep="", comment.char="#", header=FALSE)
{
    csvClasses <- xstringClasses <-
        list(sequence="DNAString", alignQuality="integer",
             nMatch="integer", position="character", strand="factor",
             refSequence=NULL, nextBestAlignQuality="integer")
    xstringNames <- "sequence"
    csvClasses[xstringNames] <- list(NULL)
    xstringClasses[!names(xstringClasses) %in% xstringNames] <-
        list(NULL)

    ## CSV portion
    lst <- .read_csv_portion(dirPath, pattern, csvClasses, ...,
                             col.names=names(csvClasses),
                             quote=quote, sep=sep,
                             comment.char=comment.char, header=header)
    idx <- regexpr(":", lst[["position"]], fixed=TRUE)
    chromosome <- substr(lst[["position"]], 1, idx-1)
    chromosome[idx==-1] <- NA
    chromosome <- factor(chromosome)
    position <- as.integer(substr(lst[["position"]], idx+1,
                                  nchar(lst[["position"]])))
    df <- data.frame(nMatch=lst$nMatch,
                     nextBestAlignQuality=lst$nextBestAlignQuality)
    meta <- data.frame(labelDescription=c(
                         "Number of matches",
                         "Next-best alignment quality score"))
    alignData <- AlignedDataFrame(df, meta)

    ## XStringSet classes
    sets <- readXStringColumns(dirPath, pattern, xstringClasses,
                               ..., sep=" \t")
    len <- length(sets[["sequence"]])
    wd <- width(sets[["sequence"]])
    q <- paste(rep(" ", max(wd)), collapse="")
    quality <- BStringSet(Views(BString(q), start=rep(1, len), end=wd))
    AlignedRead(sread=sets[["sequence"]],
                id=BStringSet(character(len)),
                quality=SFastqQuality(quality),
                chromosome=chromosome,
                position=position,
                strand=.toStrand_Solexa(lst[["strand"]]),
                alignQuality=NumericQuality(lst[["alignQuality"]]),
                alignData=alignData)
}

.readAligned_SolexaResult <-
    function(dirPath, pattern=character(0), ...,
             sep="\t", comment.char="#", quote="", header=FALSE)
{
    csvClasses <- xstringClasses <-
        list(id=NULL, sequence="DNAString", matchCode="factor",
             nExactMatch="integer", nOneMismatch="integer",
             nTwoMismatch="integer", chromosome="factor",
             position="integer", strand="factor",
             NCharacterTreatment="factor", mismatchDetailOne="character",
             mismatchDetailTwo="character")
    xstringNames <- "sequence"
    csvClasses[xstringNames] <- list(NULL)
    xstringClasses[!names(xstringClasses) %in% xstringNames] <-
        list(NULL)

    ## CSV portion
    lst <- .read_csv_portion(dirPath, pattern, csvClasses, ...,
                             col.names=names(csvClasses),
                             quote=quote, sep=sep,
                             comment.char=comment.char, header=header)
    df <- data.frame(matchCode=lst[["matchCode"]],
                     nExactMatch=lst[["nExactMatch"]],
                     nOneMismatch=lst[["nOneMismatch"]],
                     nTwoMismatch=lst[["nTwoMismatch"]],
                     NCharacterTreatment=lst[["NCharacterTreatment"]],
                     mismatchDetailOne=lst[["mismatchDetailOne"]],
                     mismatchDetailTwo=lst[["mismatchDetailTwo"]])
    meta <- data.frame(labelDescription=c(
                         "Type of match; see ?'readAligned,character-method'",
                         "Number of exact matches",
                         "Number of 1-error mismatches",
                         "Number of 2-error mismatches",
                         "Treatment of 'N'; .: NA; D: deletion; |: insertion",
                         "Mismatch error 1 detail; see ?'readAligned,character-method",
                         "Mismatch error 2 detail; see ?'readAligned,character-method"))
    alignData <- AlignedDataFrame(df, meta)
    ## XStringSet classes
    sets <- readXStringColumns(dirPath, pattern, xstringClasses,
                               ..., sep=sep)
    len <- length(sets[["sequence"]])
    wd <- width(sets[["sequence"]])
    q <- paste(rep(" ", max(wd)), collapse="")
    sfq <- BStringSet(Views(BString(q), start=rep(1, len), end=wd))
    AlignedRead(sread=sets[["sequence"]],
                quality=SFastqQuality(sfq),
                chromosome=lst[["chromosome"]],
                position=lst[["position"]],
                strand=.toStrand_Solexa(lst[["strand"]]),
                alignQuality=NumericQuality(rep(NA_integer_,
                  length(sfq))),
                alignData=alignData)
}

.SolexaExport_AlignedDataFrame <- function(data)
{
    AlignedDataFrame(data=data,
                     metadata=data.frame(
                       labelDescription=c(
                         "Analysis pipeline run",
                         "Flow cell lane",
                         "Flow cell tile",
                         "Cluster x-coordinate",
                         "Cluster y-coordinate",
                         "Read successfully passed filtering?")))
}

.readAligned_SolexaExport <-
  function(dirPath, pattern=character(0), ..., sep="\t",
           commentChar="#")
{
    files <- .file_names(dirPath, pattern)
    .Call(.read_solexa_export, files, sep, commentChar)
}

.readAligned_Maq_ADF <- function(lst) {
    df <- with(lst, data.frame(nMismatchBestHit=nMismatchBestHit,
                               mismatchQuality=mismatchQuality,
                               nExactMatch24=nExactMatch24,
                               nOneMismatch24=nOneMismatch24))
    meta <- data.frame(labelDescription=c(
                         "Number of mismatches of the best hit",
                         "Sum of mismatched base qualities of the best hit",
                         "Number of 0-mismatch hits of the first 24 bases",
                         "Number of 1-mismatch hits of the first 24 bases"))
    AlignedDataFrame(df, meta)
}


.maqmap_warning_seen <- local({
    seen <- FALSE
    function() {
        if (!seen) {
            seen <<- TRUE
            FALSE
        } else seen
    }
})

.maqmap_file_list_error <- 
    function(files)
{
    .throw(SRError("UserArgumentMismatch",
                   "%s for '%s' must match 1 file, got\n  %s",
                   "'dirPath', 'pattern'", "MAQMap",
                   paste(files, collapse="\n  ")))
}

.readAligned_MaqMapOld <-
    function(dirPath, pattern=character(0), records=-1L, ...)
{
    files <- .file_names(dirPath, pattern)
    if (length(files) > 1)
        .maqmap_file_list_error(files)
    lst <- .Call(.read_maq_map, files, as.integer(records), FALSE)
    AlignedRead(sread=lst[["readSequence"]],
                id=lst[["readId"]],
                quality=FastqQuality(lst[["fastqScores"]]),
                chromosome=lst[["chromosome"]],
                position=lst[["position"]],
                strand=lst[["strand"]],
                alignQuality=IntegerQuality(lst[["alignQuality"]]),
                alignData=.readAligned_Maq_ADF(lst))
}
        
.readAligned_MaqMap <-
    function(dirPath, pattern=character(0), records=-1L, ...)
{
    files <- .file_names(dirPath, pattern)
    if (length(files) > 1)
        .maqmap_file_list_error(files)
    if (!.maqmap_warning_seen()) {
        msg <- paste("API change: The type 'MAQMap' now reads map files produced",
                     "with at least version 0.7.0 of Maq. Before version 0.99.3 of the",
                     "ShortRead package, this type was for reading map files produced by",
                     "Maq up to version 0.6.x. If you still use the old Maq version, change",
                     "the type argument to 'MAQMapShort'. [This warning will disappear in a",
                     "future version of ShortRead.]")
        warning( paste(strwrap(msg, indent=2, exdent=2), collapse="\n") )
    }
    lst <- .Call(.read_maq_map, files, as.integer(records), TRUE)
    AlignedRead(sread=lst[["readSequence"]],
                id=lst[["readId"]],
                quality=FastqQuality(lst[["fastqScores"]]),
                chromosome=lst[["chromosome"]],
                position=lst[["position"]],
                strand=lst[["strand"]],
                alignQuality=IntegerQuality(lst[["alignQuality"]]),
                alignData=.readAligned_Maq_ADF(lst))
}

.readAligned_MaqMapview <-
    function(dirPath, pattern=character(0), ..., sep="\t", header=FALSE,
             quote="")
{
    colClasses <-
        list(NULL, chromosome="factor", position="integer",
             strand="factor", NULL, NULL, alignQuality="integer", NULL,
             NULL, nMismatchBestHit="integer",
             mismatchQuality="integer", nExactMatch24="integer",
             nOneMismatch24="integer", NULL, NULL, NULL)
    ## CSV portion
    csv <- .read_csv_portion(dirPath, pattern, colClasses, sep=sep,
                             header=header, quote=quote, ...)
    ## XStringSet components
    colClasses <- list("BString", NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, "DNAString", "BString")
    sets <- readXStringColumns(dirPath, pattern,
                               colClasses, sep=sep, header=header)
    AlignedRead(sread=sets[[2]], id=sets[[1]],
                quality=FastqQuality(sets[[3]]),
                chromosome=csv[["chromosome"]],
                position=csv[["position"]],
                strand=factor(csv[["strand"]], levels=.STRAND_LEVELS),
                alignQuality=IntegerQuality(csv[["alignQuality"]]),
                alignData=.readAligned_Maq_ADF(csv))
}

.Bowtie_AlignedDataFrame <- function(mismatch)
{
    df <- data.frame(mismatch=mismatch, stringsAsFactors=FALSE)
    meta <- data.frame(labelDescription=c(
                         "Comma-separated mismatch positions"))
    AlignedDataFrame(df, meta)
}

.readAligned_Bowtie <-
  function(dirPath, pattern=character(0), ...,
           qualityType="SFastqQuality", sep="\t", commentChar="#")
{
    files <- .file_names(dirPath, pattern)
    .Call(.read_bowtie, files, qualityType, sep, commentChar)
}

.SOAP_AlignedDataFrame <-
    function(nEquallyBestHits, pairedEnd, alignedLength,
             typeOfHit, hitDetail)
{
    df <- data.frame(nEquallyBestHits=nEquallyBestHits,
                     pairedEnd=factor(pairedEnd),
                     alignedLength=alignedLength,
                     typeOfHit=typeOfHit,
                     hitDetail=hitDetail,
                     stringsAsFactors=FALSE)
    meta <- data.frame(labelDescription=c(
                         "Number of equally-best hits",
                         "Paired end, a or b",
                         "Length of read used in alignment",
                         "Integer indicator of match type; 0: exact; 1-100: mismatch; 100+n: n-base insertion; 200+n: n-base deletion", 
                         "Detailed description of match"))
    AlignedDataFrame(df, meta)
}

.readAligned_SOAP <-
    function(dirPath, pattern=character(0), ...,
             qualityType="SFastqQuality", sep="\t", commentChar="#")
{
    files <- .file_names(dirPath, pattern)
    .Call(.read_soap, files, qualityType, sep, commentChar)
}

.readAligned_character <-
    function(dirPath, pattern=character(0),
             type=c(
               "SolexaExport", "SolexaAlign",
               "SolexaPrealign", "SolexaRealign",
               "SolexaResult",
               "MAQMap", "MAQMapShort", "MAQMapview",
               "Bowtie", "SOAP"),
             ..., filter=srFilter())
{
    if (missing(type))
        .arg_missing_err("type", "readAligned,character-method",
                       "help(\"readAligned,character-method\")")
    if (!is.character(type) || length(type) != 1)
        .arg_mismatch_type_err("type", "character(1)")
    if (!missing(filter))
        .check_type_and_length(filter, "SRFilter", NA)
    vals <- eval(formals(sys.function())$type)
    if (!type %in% vals)
        .arg_mismatch_value_err("type", type, vals)
    aln <- 
        tryCatch({
            switch(type,
                   SolexaExport=.readAligned_SolexaExport(dirPath,
                     pattern=pattern, ...),
                   SolexaPrealign=,
                   SolexaAlign=,
                   SolexaRealign=.readAligned_SolexaAlign(dirPath,
                     pattern=pattern, ...),
                   SolexaResult=.readAligned_SolexaResult(dirPath,
                     pattern=pattern, ...),
                   MAQMap=.readAligned_MaqMap(dirPath, pattern, ...),
                   MAQMapShort=.readAligned_MaqMapOld(dirPath, pattern, ...),
                   MAQMapview=.readAligned_MaqMapview(
                     dirPath, pattern=pattern, ...),
                   Bowtie=.readAligned_Bowtie(dirPath, pattern, ...),
                   SOAP=.readAligned_SOAP(dirPath, pattern, ...))
        }, error=function(err) {
            if (is(err, "SRError")) stop(err)
            else {
                pat <- paste(pattern, collapse=" ")
                txt <- paste("'%s' failed to parse files",
                             "dirPath: '%s'",
                             "pattern: '%s'",
                             "type: '%s'",
                             "error: %s", sep="\n  ")
                msg <- sprintf(txt, "readAligned",
                               paste(dirPath, collapse="'\n    '"),
                               paste(pat, collapse="'\n    '"),
                               type, conditionMessage(err))
                .throw(SRError("Input/Output", msg))
            }
        })
    if (!missing(filter))
        aln <- aln[filter(aln)]
    aln
}

setMethod(readAligned, "character", .readAligned_character)
