.AlignedRead_validity <- function(object) {
    msg <- NULL
    len <- length(sread(object))
    slts <- c("chromosome", "position", "strand", "alignQuality")
    olen <- sapply(slts, function(elt) {
        length(do.call(elt, list(object)))
    })
    if (!all(olen==len)) {
        bad <- olen!=len
        msg <- c(msg,
                 sprintf("length mismatch: expected %d, found:\n  %s",
                         len, paste(slts[bad], olen[bad], sep="=",
                                    collapse=", ")))
    }
    if (is.null(msg)) TRUE else msg
}

setMethod(".srValidity", "AlignedRead", .AlignedRead_validity)

AlignedRead <- function(sread, id, quality,
                        chromosome, position, strand,
                        alignQuality,
                        alignData=AlignedDataFrame(nrow=length(sread))) {
    new("AlignedRead",
        sread=sread, id=id, quality=quality,
        chromosome=chromosome, position=position, strand=strand,
        alignQuality=alignQuality, alignData=alignData)
}

.read_csv_portion <- function(dirPath, pattern, colClasses, sep, header) {
    ## visit several files, then collapse
    files <- list.files(dirPath, pattern, full.names=TRUE)
    if (length(files)==0)
        .throw(SRError("Input/Output",
                       "no input files found\n  dirPath: %s\n  pattern: %s\n",
                       dirPath, pattern))
    lsts <- lapply(files, read.csv,
                   sep=sep, header=header, colClasses=colClasses)
    cclasses <- colClasses[!sapply(colClasses, is.null)]
    lst <- lapply(seq_along(names(cclasses)),
                  function(idx) unlist(lapply(lsts, "[[", idx)))
    names(lst) <- names(cclasses)
    lst
}

.readAligned_SolexaExport <- function(dirPath, pattern=character(0),
                                      sep="\t", header=FALSE) {
    ## NULL are currently ignored, usually from paired-end reads
    csvClasses <- xstringClasses <-
        list(machine=NULL, run="integer", lane="integer",
             tile="integer", x="integer", y="integer",
             indexString=NULL, pairedReadNumber=NULL,
             sequence="DNAString", quality="BString",
             chromosome="factor", contig=NULL, position="integer",
             strand="factor", descriptor=NULL, alignQuality="integer",
             pairedScore=NULL, partnerCzome=NULL, partnerContig=NULL,
             partnerOffset=NULL, partnerStrand=NULL,
             filtering="factor")

    xstringNames <- c("sequence", "quality")
    csvClasses[xstringNames] <- list(NULL, NULL)
    xstringClasses[!names(xstringClasses) %in% xstringNames] <-
        list(NULL)

    ## CSV portion
    lst <- .read_csv_portion(dirPath, pattern, csvClasses, sep, header)
    df <- with(lst, data.frame(run=run, lane=lane, tile=tile, x=x,
                               y=y, filtering=filtering))
    meta <- data.frame(labelDescription=c(
                         "Analysis pipeline run",
                         "Flow cell lane",
                         "Flow cell tile",
                         "Cluster x-coordinate",
                         "Cluster y-coordinate",
                         "Read successfully passed filtering?"))
    alignData <- AlignedDataFrame(df, meta)

    ## XStringSet classes
    sets <- readXStringColumns(dirPath, pattern, xstringClasses,
                               sep=sep, header=header)

    AlignedRead(sread=sets[["sequence"]],
                id=BStringSet(character(length(sets[["sequence"]]))),
                quality=SFastqQuality(sets[["quality"]]),
                chromosome=lst[["chromosome"]],
                position=lst[["position"]],
                strand=lst[["strand"]],
                alignQuality=NumericQuality(lst[["alignQuality"]]),
                alignData=alignData)
}

.readAligned_MaqMapview <- function(dirPath, pattern=character(0), 
                                    sep="\t", header=FALSE) {
    colClasses <-
        list(NULL, chromosome="factor", position="integer",
             strand="factor", NULL, NULL, quality="numeric", NULL,
             NULL, nMismatchBestHit="integer",
             mismatchQuality="numeric", nExactMatch24="integer",
             nOneMismatch24="integer", NULL, NULL, NULL)

    ## CSV portion
    lst <- .read_csv_portion(dirPath, pattern, colClasses, sep,
                             header)
    df <- with(lst, data.frame(nMismatchBestHit=nMismatchBestHit,
                               mismatchQuality=mismatchQuality,
                               nExactMatch24=nExactMatch24,
                               nOneMismatch24=nOneMismatch24))
    meta <- data.frame(labelDescription=c(
                         "Number of mismatches of the best hit",
                         "Sum of mismatched base qualities of the best hit",
                         "Number of 0-mismatch hits of the first 24 bases",
                         "Number of 1-mismatch hits of the first 24 bases"))
    alignData <- AlignedDataFrame(df, meta)

    ## XStringSet components
    colClasses <- list("BString", NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, "DNAString", "BString")
    sets <- readXStringColumns(dirPath, pattern,
                               colClasses, sep=sep, header=header)

    AlignedRead(sread=sets[[2]], id=sets[[1]],
                quality=FastqQuality(sets[[3]]),
                chromosome=lst[["chromosome"]],
                position=lst[["position"]],
                strand=lst[["strand"]],
                alignQuality=NumericQuality(lst[["quality"]]),
                alignData=alignData)
}

.readAligned_character<- function(dirPath, pattern=character(0),
                                  type=c("SolexaExport",
                                    "MAQMapview"), ...) {
  if (!is.character(type) || length(type) != 1)
    .throw(SRError("UserArgumentMismatch",
                   "'%s' must be '%s'",
                   "type", "character(1)"))
  switch(type,
         SolexaExport=.readAligned_SolexaExport(dirPath,
           pattern=pattern, ...),
         MAQMapview=.readAligned_MaqMapview(dirPath, pattern=pattern,
           ...),
         .throw(SRError("UserArgumentMismatch",
                        "'%s' unknown; value was '%s'",
                        "type", type)))
}

setMethod("readAligned", "character", .readAligned_character)

.make_getter(c("chromosome", "position", "strand", "alignQuality",
               "alignData"))

## subset

setMethod("[", c("AlignedRead", "missing", "missing"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("AlignedRead", "missing", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("AlignedRead", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

.AlignedRead_subset <- function(x, i, j, ..., drop=TRUE) {
    if (nargs() != 2) .subset_err()
    initialize(x, sread=sread(x)[i], id=id(x)[i], quality=quality(x)[i],
               chromosome=chromosome(x)[i], position=position(x)[i],
               strand=strand(x)[i], alignQuality=alignQuality(x)[i],
               alignData=alignData(x)[i,])
}

setMethod("[", c("AlignedRead", "ANY", "missing"), .AlignedRead_subset)

## show

setMethod("show", "AlignedRead", function(object) {
    callNextMethod()
    cat("chromosome:", selectSome(chromosome(object)), "\n")
    cat("position:", selectSome(position(object)), "\n")
    cat("strand:", selectSome(strand(object)), "\n")
    cat("alignQuality:", class(alignQuality(object)), "\n")
    cat("alignData varLabels:",
        selectSome(varLabels(alignData(object))), "\n")
})

setMethod("detail", "AlignedRead", function(object, ...) {
    callNextMethod()
    cat("\nchromosome:", selectSome(chromosome(object)), "\n")
    cat("position:", selectSome(position(object)), "\n")
    cat("strand:", selectSome(strand(object)), "\n")
    cat("alignQuality:\n")
    detail(alignQuality(object))
    cat("\nalignData:\n")
    show(alignData(object))
})
