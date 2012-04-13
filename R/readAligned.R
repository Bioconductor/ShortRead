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
             quote="", sep="", comment.char="#", header=FALSE,
             Lpattern="", Rpattern="")
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
    lbls <- c(run="Analysis pipeline run",
              lane="Flow cell lane",
              tile="Flow cell tile",
              x="Cluster x-coordinate",
              y="Cluster y-coordinate",
              filtering="Read successfully passed filtering?",
              contig="Contig",
              multiplexIndex="Multiplex index",
              pairedReadNumber="Paired read number")[names(data)]
    AlignedDataFrame(data=data, metadata=data.frame(labelDescription=lbls))
}

.readAligned_SolexaExport <-
  function(dirPath, pattern=character(0), ..., withAll=FALSE,
           withId=withAll, withMultiplexIndex=withAll,
           withPairedReadNumber=withAll,
           sep="\t", commentChar="#")
{
    files <- .file_names(dirPath, pattern)
    .Call(.read_solexa_export, files, sep, commentChar,
          c(withId, withMultiplexIndex, withPairedReadNumber))
}

.readAligned_Maq_ADF <- function(lst) {
    df <- data.frame(nMismatchBestHit=lst$nMismatchBestHit,
                     mismatchQuality=lst$mismatchQuality,
                     nExactMatch24=lst$nExactMatch24,
                     nOneMismatch24=lst$nOneMismatch24)
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
    function(files, type)
{
    .throw(SRError("UserArgumentMismatch",
                   "%s for '%s' must match 1 file, got\n  %s",
                   "'dirPath', 'pattern'", type,
                   paste(files, collapse="\n  ")))
}

.readAligned_MaqMapOld <-
    function(dirPath, pattern=character(0), records=-1L, ...)
{
    files <- .file_names(dirPath, pattern)
    if (length(files) > 1)
        .maqmap_file_list_error(files, "MAQMapShort")
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
        .maqmap_file_list_error(files, "MAQMap")
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

.Bowtie_AlignedDataFrame <- function(similar, mismatch)
{
    df <- data.frame(similar=similar, mismatch=mismatch,
                     stringsAsFactors=FALSE)
    meta <- data.frame(labelDescription=c(
                         "if Bowtie >= 0.9.9.3 (May 12, 2009)?: number of alignments aligning to the same reference characters; else 'Reserved'",
                         "Comma-separated mismatch positions"))
    AlignedDataFrame(df, meta)
}

.readAligned_Bowtie <-
  function(dirPath, pattern=character(0), ...,
           qualityType=c("FastqQuality", "SFastqQuality"),
           sep="\t", commentChar="#")
{
    tryCatch(qualityType <- match.arg(qualityType),
             error=function(err) {
                 .throw(SRError("UserArgumentMismatch",
                                conditionMessage(err)))
             })
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

.readAligned_bamWhat <- function(withQname=TRUE)
{
    c(if (withQname) "qname" else NULL,
      "flag", "rname", "strand", "pos", "mapq", "seq", "qual")
}

.readAligned_bam <-
    function(dirPath, pattern=character(0), ...,
             param=ScanBamParam(simpleCigar=TRUE,
               reverseComplement=TRUE, what=.readAligned_bamWhat()))
{
    files <-
        if (!all(grepl("^(ftp|http)://", dirPath)))
            .file_names(dirPath, pattern)
        else {
            if (length(dirPath) != 1 || length(pattern) != 0) {
                msg <- paste("ftp:// and http:// support requires",
                             "'dirPath' as character(1),",
                             "'pattern' as character(0)", collapse="")
                .throw(SRError("UserArgumentMismatch", msg))
            }
            dirPath
        }

    if (!is(param, "ScanBamParam")) {
        msg <-  "'param' must be a ScanBamParam object."
        .throw(SRError("UserArgumentMismatch",msg))
    }
    ## FIXME: currently we only deal with cigars without indels
    if (bamSimpleCigar(param) != TRUE) {
        msg <- paste("using 'TRUE' for 'bamSimpleCigar(param)'",
                     "(skipping reads with I, D, H, S, or P in 'cigar')")
        .throw(SRWarn("UserArgumentMismatch", msg))
        bamSimpleCigar(param) <- TRUE
    }
    if (bamReverseComplement(param) != TRUE) {
        msg <- "using 'TRUE' for 'bamReverseComplement(param)'"
        .throw(SRWarn("UserArgumentMismatch", msg))
        bamReverseComplement(param) <- TRUE
    }
    what <- .readAligned_bamWhat("qname" %in% bamWhat(param))
    if (!(length(bamWhat(param)) && all(what %in% bamWhat(param)))) {
        msg <- sprintf("using at least '%s' for 'bamWhat(param)'",
                       paste(what, collapse="' '"))
        .throw(SRWarn("UserArgumentMismatch", msg))
        bamWhat(param) <- union(bamWhat(param), what)
    }

    ## handle multiple files and params
    result <- mapply(scanBam, files, MoreArgs=list(param=param),
                     ..., SIMPLIFY=FALSE, USE.NAMES=FALSE)
    ulist <- function(X, ..., recursive=TRUE)
        unlist(lapply(X, lapply, "[[", ...),
               recursive=recursive, use.names=FALSE)
    cxslist <- function(X, ...)
        do.call(c, ulist(X, ...))
    chromosome <- local({
        X <- ulist(result, "rname", recursive=FALSE)
        values <- do.call(c, lapply(X, as.character))
        factor(values, levels=unique(unlist(lapply(X, levels))))
    })
    strand <- local({
        X <- ulist(result, "strand", recursive=FALSE)
        values <- do.call(c, lapply(X, as.character))
        strand(values)
    })
    id <- local({
        X <- ulist(result, "qname")
        if (is.null(X)) BStringSet(character(length(strand)))
        else BStringSet(X)
    })

    AlignedRead(sread=cxslist(result, "seq"), id=id,
                quality=FastqQuality(as(cxslist(result, "qual"),
                  "BStringSet")),
                chromosome=chromosome, strand=strand,
                position=ulist(result, "pos"),
                alignQuality=NumericQuality(ulist(result, "mapq")),
                alignData=AlignedDataFrame(
                  data=data.frame(flag=ulist(result, "flag")),
                  metadata=data.frame(
                    labelDescription=c("Type of read; see ?scanBam"))))
}

.readAligned_character <-
    function(dirPath, pattern=character(0),
             type=c(
               "SolexaExport", "SolexaAlign",
               "SolexaPrealign", "SolexaRealign",
               "SolexaResult",
               "MAQMap", "MAQMapShort", "MAQMapview",
               "Bowtie", "SOAP", "BAM"),
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
                   SOAP=.readAligned_SOAP(dirPath, pattern, ...),
                   BAM=.readAligned_bam(dirPath, pattern, ...))
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
