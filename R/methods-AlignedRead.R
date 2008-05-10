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

.readMaqMapview <- function(dirPath, pattern=character(0), 
                           sep="\t", header=FALSE) {
    colClasses <- list(NULL, "factor", "integer", "factor", NULL,
                       NULL, "numeric", NULL, NULL, "integer",
                       "numeric", "integer", "integer", NULL, NULL,
                       NULL)
    cNames <- c("chromosome", "position", "strand", "quality",
                "nMismatchBestHit", "mismatchQuality", "nExactMatch24",
                "nOneMismatch24")
    names(colClasses)[!sapply(colClasses, is.null)] <- cNames
    ## visit several files, then collapse
    files <- list.files(dirPath, pattern, full.names=TRUE)
    lsts <- lapply(files, read.csv,
                   sep=sep, header=header, colClasses=colClasses)
    cclasses <- colClasses[!sapply(colClasses, is.null)]
    lst <- lapply(seq_along(cNames),
                  function(idx) unlist(lapply(lsts, "[[", idx)))
    names(lst) <- cNames
    ## alignedData
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

.readAligned_character<- function(dirPath,
                                  pattern=character(0),
                                  type="MAQMapview", ...) {
  if (!is.character(type) || length(type) != 1)
    .throw(SRError("UserArgumentMismatch",
                   "'%s' must be '%s'",
                   "type", "character(1)"))
  switch(type,
         MAQMapview=.readMaqMapview(dirPath, pattern=pattern, ...),
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
    cat("chromosome:", selectSome(levels(chromosome(object))), "\n")
    cat("position:", selectSome(position(object)), "\n")
    cat("strand:", selectSome(as.character(strand(object))), "\n")
    cat("alignQuality:", class(alignQuality(object)), "\n")
    cat("alignData varLabels:",
        selectSome(varLabels(alignData(object))), "\n")
})

setMethod("detail", "AlignedRead", function(object, ...) {
    callNextMethod()
    cat("\nchromosome:", selectSome(levels(chromosome(object))), "\n")
    cat("position:", selectSome(position(object)), "\n")
    cat("strand:", selectSome(as.character(strand(object))), "\n")
    cat("alignQuality:\n")
    detail(alignQuality(object))
    cat("\nalignData:\n")
    show(alignData(object))
})
