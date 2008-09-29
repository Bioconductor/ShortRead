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

.make_getter(c("chromosome", "position", "alignQuality",
               "alignData"))

setMethod("strand", "AlignedRead", function(object, ...) 
{
    slot(object, "strand")
})

## coerce

setAs("PairwiseAlignment", "AlignedRead",
      function(from, to) {
        pat <- pattern(from)
        quality <- character()
        if (is(pat, "QualityAlignedXStringSet"))
          quality <- quality(pat)
        new("AlignedRead", sread = unaligned(pat), id = names(pat),
            quality = FastqQuality(quality),
            position = start(Views(pat)),
            alignQuality = IntegerQuality(score(from)))
      })

## subset

setMethod("[", c("AlignedRead", "missing", "missing"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("AlignedRead", "missing", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("AlignedRead", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

.AlignedRead_subset <- function(x, i, j, ..., drop=TRUE) {
    if (nargs() != 2) .subset_err()
    initialize(x, sread=sread(x)[i], id=id(x)[i],
               quality=quality(x)[i], chromosome=chromosome(x)[i],
               position=position(x)[i], strand=strand(x)[i],
               alignQuality=alignQuality(x)[i],
               alignData=alignData(x)[i,]) }

setMethod("[", c("AlignedRead", "ANY", "missing"), .AlignedRead_subset)

## manip

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

## summary

## perhaps summary statistics like ShortReadQ except broken down by chromosome,
## strand, and their combination
