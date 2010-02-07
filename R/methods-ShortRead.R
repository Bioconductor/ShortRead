.ShortRead_validity <- function(object) {
    msg <- NULL
    if (length(sread(object)) != length(id(object)))
        msg <- c(msg,
                 sprintf("sread() and id() length mismatch: %d, %d",
                         length(sread(object)), length(id(object))))
    if (is.null(msg)) TRUE else msg
}

setMethod(.srValidity, "ShortRead", .ShortRead_validity)

.make_getter("sread")

setMethod(id, "ShortRead",
          function(object, ...) slot(object, "id"))

setMethod(ShortRead, c("DNAStringSet", "BStringSet"),
          function(sread, id, ...)
{
    new("ShortRead", sread=sread, id=id, ...)
})

setMethod(ShortRead, c("DNAStringSet", "missing"),
          function(sread, id, ...)
{
    new("ShortRead", sread=sread, id=BStringSet(rep("", length(sread))), ...)
})

setMethod(ShortRead, c("missing", "missing"),
          function(sread, id, ...)
{
    new("ShortRead")
})

setMethod(length, "ShortRead", function(x) length(sread(x)))

setMethod(width, "ShortRead", function(x) width(sread(x)))

## coerce

setMethod(pairwiseAlignment, "ShortRead",
          function(pattern, subject, ...)
          {
            pairwiseAlignment(sread(pattern), subject, ...)
          })

## import

setMethod(readFasta, "character", function(dirPath, pattern=character(),
                                             sample = 1, ...) {
  src <- .file_names(dirPath, pattern)[sample]
  FASTAlist <- lapply(src, readFASTA, strip.desc = TRUE)
  FASTArecs <- do.call(c, FASTAlist)
  strings <- FASTArecordsToCharacter(FASTArecs)
  new("ShortRead", ...,
      sread=DNAStringSet(strings, use.names=FALSE),
      id=BStringSet(names(strings)))
})

## subset

setMethod("[", c("ShortRead", "missing", "missing"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("ShortRead", "missing", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("ShortRead", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

.ShortRead_subset <- function(x, i, j, ..., drop=TRUE) {
    if (nargs() != 2) .subset_err()
    initialize(x, sread=sread(x)[i], id=id(x)[i])
}

setMethod("[", c(x="ShortRead", i="ANY", j="missing"),
          .ShortRead_subset)

setMethod(append, c("ShortRead", "ShortRead", "missing"),
    function(x, values, after=length(x)) 
{
    initialize(x, id=append(id(x), id(values)),
               sread=append(sread(x), sread(values)))
})

## manip

.abc_ShortRead <- function(stringSet, alphabet, ...) {
    if (!missing(alphabet))
        .throw(SRWarn("UserArgumentMismatch", "'alphabet' ignored"))
    alphabetByCycle(sread(stringSet), ...)
}

setMethod(alphabetByCycle, "ShortRead", .abc_ShortRead)

setMethod(clean, "ShortRead", function(object, ...) {
    alf <- alphabetFrequency(sread(object), baseOnly=TRUE)
    object[alf[,'other'] == 0]
})

setMethod(dustyScore, "ShortRead",
          function(x, ...) callGeneric(sread(x), ...))

setMethod(srorder, "ShortRead", .forward_x)

setMethod(srrank, "ShortRead", .forward_x)

setMethod(srsort, "ShortRead", function(x, ...) {
    x[srorder(x, ...)]
})

setMethod(srduplicated, "ShortRead", .forward_x)

setMethod(tables, "ShortRead", function(x, n=50, ...) {
    callGeneric(sread(x), n=n, ...)
})

.srdistance_ShortRead_ANY <- function(pattern, subject, ...)
{
    callGeneric(sread(pattern), subject, ...)
}

setMethod(srdistance, c("ShortRead", "ANY"),
          .srdistance_ShortRead_ANY)

setMethod(narrow, "ShortRead",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
{
    initialize(x,
               sread=narrow(sread(x), start, end, width, use.names))
})

setMethod(compact, "ShortRead",
    function(x, ...)
{
    initialize(x, id=callGeneric(id(x)), sread=callGeneric(sread(x)))
})

setMethod(trimLRPatterns, c(subject="ShortRead"),
    function (Lpattern = "", Rpattern = "", subject, max.Lmismatch =
              0, max.Rmismatch = 0, with.Lindels = FALSE, with.Rindels
              = FALSE, Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
{
    ranges <-
        callGeneric(Lpattern, Rpattern, sread(subject), max.Lmismatch,
                    max.Rmismatch, with.Lindels, with.Rindels, Lfixed,
                    Rfixed, ranges=TRUE)
    narrow(subject, start(ranges), end(ranges))
})

## show

setMethod(show, "ShortRead", function(object) {
    callNextMethod()
    wd <- sort(unique(width(object)))
    if (length(wd)>2) wd <- paste(range(wd), collapse="..")
    cat("length:", length(object), "reads; width:", wd, "cycles\n")
})

setMethod(detail, "ShortRead", function(object, ...) {
    cat("class: ", class(object), "\n")
    cat("\nsread:\n")
    show(sread(object))
    cat("\nid:\n")
    show(id(object))
})

## summary

## perhaps a 'summary' method with statistics on reads for each sample
