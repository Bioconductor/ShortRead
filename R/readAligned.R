.read_csv_portion <- function(dirPath, pattern, colClasses, ...) {
    ## visit several files, then collapse
    files <- .file_names(dirPath, pattern)
    lsts <- lapply(files, read.csv,
                   ..., colClasses=colClasses)
    cclasses <- colClasses[!sapply(colClasses, is.null)]
    lst <- lapply(seq_along(names(cclasses)),
                  function(idx) unlist(lapply(lsts, "[[", idx)))
    names(lst) <- names(cclasses)
    lst
}

.readAligned_SolexaAlign <-
    function(dirPath, pattern=character(0), ...,
             sep=" ", comment.char="#", header=FALSE)
{
    csvClasses <- xstringClasses <-
        list(sequence="DNAString", alignQuality="integer",
             nMatch="integer", position="integer", strand="factor",
             refSequence=NULL, nextBestAlignQuality="integer")
    xstringNames <- "sequence"
    csvClasses[xstringNames] <- list(NULL)
    xstringClasses[!names(xstringClasses) %in% xstringNames] <-
        list(NULL)

    ## CSV portion
    lst <- .read_csv_portion(dirPath, pattern, csvClasses, ...,
                             sep=sep, comment.char=comment.char,
                             header=header)
    df <- data.frame(nMatch=lst$nMatch,
                     nextBestAlignQuality=lst$nextBestAlignQuality)
    meta <- data.frame(labelDescription=c(
                         "Number of matches",
                         "Next-best alignment quality score"))
    alignData <- AlignedDataFrame(df, meta)

    ## XStringSet classes
    sets <- readXStringColumns(dirPath, pattern, xstringClasses,
                               ..., sep=sep)
    len <- length(sets[["sequence"]])
    wd <- width(sets[["sequence"]])
    q <- paste(rep(" ", max(wd)), collapse="")
    quality <- BStringSet(Views(BString(q), start=rep(1, len), end=wd))
    AlignedRead(sread=sets[["sequence"]],
                id=BStringSet(character(len)),
                quality=SFastqQuality(quality),
                chromosome=factor(rep(NA, len)),
                position=lst[["position"]],
                strand=factor(lst[["strand"]], levels=.STRAND_LEVELS),
                alignQuality=NumericQuality(lst[["alignQuality"]]),
                alignData=alignData)
}

.SolexaExport_AlignedDataFrame <-
    function(run, lane, tile, x, y, filtering)
{
    df <- data.frame(run=factor(run),
                     lane=lane, tile=tile, x=x, y=y,
                     filtering=filtering)
    meta <- data.frame(labelDescription=c(
                         "Analysis pipeline run",
                         "Flow cell lane",
                         "Flow cell tile",
                         "Cluster x-coordinate",
                         "Cluster y-coordinate",
                         "Read successfully passed filtering?"))
    AlignedDataFrame(df, meta)
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

.readAligned_MaqMapOld <-
    function(dirPath, pattern=character(0), records=-1L, ...)
{
    files <- .file_names(dirPath, pattern)
    if (length(files) > 1)
        .arg_mismatch_type_err("dirPath', 'pattern", "character(1)")
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
        .arg_mismatch_type_err("dirPath', 'pattern", "character(1)")
    warning( paste( "API change: The type 'MAQMap' now reads map files produced",
       "with at least version 0.7.0 of Maq. Before version 0.99.3 of the",
       "ShortRead package, this type was for reading map files produced by",
       "Maq up to version 0.6.x. If you still use the old Maq version, change",
       "the type argument to 'MAQMapShort'. [This warning will disappear in a",
       "future version of ShortRead.]" ) )
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

.readAligned_character <-
    function(dirPath, pattern=character(0),
             type=c(
               "SolexaExport", "SolexaAlign",
               "SolexaPrealign", "SolexaRealign",
               "MAQMap", "MAQMapShort", "MAQMapview"),
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
                   MAQMap=.readAligned_MaqMap(dirPath, pattern, ...),
                   MAQMapShort=.readAligned_MaqMapOld(dirPath, pattern, ...),
                   MAQMapview=.readAligned_MaqMapview(
                     dirPath, pattern=pattern, ...))
        }, error=function(err) {
            if (is(err, "SRError")) stop(err)
            else {
                pat <- paste(pattern, collapse=" ")
                txt <- paste("'%s' failed to parse files",
                             "dirPath: '%s'",
                             "pattern: '%s'",
                             "type: '%s'",
                             "error: %s", sep="\n  ")
                msg <- sprintf(txt, "readAligned", dirPath, pat,
                               type, conditionMessage(err))
                .throw(SRError("Input/Output", msg))
            }
        })
    if (!missing(filter))
        aln <- aln[filter(aln)]
    aln
}

setMethod("readAligned", "character", .readAligned_character)
