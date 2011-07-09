.abc_BStringSet <- function(stringSet, alphabet, ...) {
    if (missing(alphabet))
        alphabet <- sapply(as.raw(1:255), rawToChar)
    callNextMethod(stringSet, alphabet=alphabet)
}

setMethod(clean, "DNAStringSet", function(object, ...) {
    object[alphabetFrequency(object, baseOnly=TRUE)[,'other']==0]
})

setMethod(dustyScore, "DNAStringSet",
          function(x, batchSize=NA, ...)
{
    doDusty <- function(tripletPDict, x) {
        tnf <- t(vcountPDict(tripletPDict, x)) - 1L
        tnf[tnf < 0] <- 0L
        rowSums(tnf * tnf)    
    }
    triplets <- DNAStringSet(mkAllStrings(c("A", "C", "G", "T"), 3))
    tripletPDict <- PDict(triplets)

    if (is.na(batchSize) || length(x) <= batchSize)
        return(doDusty(tripletPDict, x))

    n <- as.integer(1L + length(x) / batchSize)
    i <- seq_len(length(x))
    i <- split(i, cut(i, n, labels=FALSE))

    unlist(unname(srapply(i, function(idx, tripletPDict, x, ...) {
        doDusty(tripletPDict, x[idx])
    }, tripletPDict, x)))
})

setMethod(alphabetByCycle, "BStringSet", .abc_BStringSet)

setMethod(srorder, "XStringSet", function(x, ...) {
    if (length(list(...))!=0)
        .throw(SRError("UserArgumentMismatch",
                       "argument '%s' not supported",
                       names(list(...))))
    .Call(.alphabet_order, x)
})

setMethod(srrank, "XStringSet", function(x, ...) {
    if (length(list(...))!=0)
        .throw(SRError("UserArgumentMismatch",
                       "argument '%s' not supported",
                       names(list(...))))
    .Call(.alphabet_rank, x)
})

setMethod(srsort, "XStringSet", function(x, ...) x[srorder(x, ...)])

setMethod(srduplicated, "XStringSet", function(x, ...) {
    if (length(list(...))!=0)
        .throw(SRError("UserArgumentMismatch",
                       "argument '%s' not supported",
                       names(list(...))))
    .Call(.alphabet_duplicated, x)
})

## srdistance

.srdistance <-
    function(pattern, subject, distanceFunc,..., verbose=FALSE)
{
    if (verbose)
        cat(".srdistance", as.character(subject), "\n")
    substitutionMatrix <- distanceFunc(pattern, subject)
    -pairwiseAlignment(pattern, subject,
                       substitutionMatrix=substitutionMatrix,
                       gapOpening=0, gapExtension=-1, scoreOnly=TRUE)
}

.srdistanceDNA <- function(pattern, subject)
{
    m <- matrix(c(1,0,0,0,.5,.5,.5,.0,.0,.0,.3,.3,.3,.0,.25,.25,.25,
                  0,1,0,0,.5,.0,.0,.5,.5,.0,.3,.3,.0,.3,.25,.25,.25,
                  0,0,1,0,.0,.5,.0,.5,.0,.5,.3,.0,.3,.3,.25,.25,.25,
                  0,0,0,1,.0,.0,.5,.0,.5,.5,.0,.3,.3,.3,.25,.25,.25),
                nrow=4, byrow=TRUE,
                dimnames=list(DNA_ALPHABET[1:4], DNA_ALPHABET))
    patternAlf <- alphabetFrequency(pattern, collapse=TRUE)
    subjectAlf <- alphabetFrequency(subject)
    alf <- unique(c(names(patternAlf)[patternAlf!=0],
                    names(subjectAlf)[subjectAlf!=0]))
    m <- m[, alf]
    -(1 - t(m) %*% m)
}

.srdistance_DNAStringSet_character <- function(pattern, subject, ...)
{
    strings <- lapply(subject, DNAString)
    res <- srapply(strings, .srdistance, pattern=pattern,
                   distanceFunc=.srdistanceDNA, ...)
    if (length(res) == length(subject))
        names(res) <- subject
    res
}

setMethod(srdistance, c("DNAStringSet", "character"),
          .srdistance_DNAStringSet_character)

.srdistance_DNAStringSet_DNAString <- function(pattern, subject, ...)
{
    res <- list(.srdistance(pattern, subject, .srdistanceDNA, ...))
    names(res) <- as.character(subject)
    res
}

setMethod(srdistance, c("DNAStringSet", "DNAString"),
          .srdistance_DNAStringSet_DNAString)

.srdistance_DNAStringSet_DNAStringSet <- function(pattern, subject,
                                                  ...)
{
    callGeneric(pattern, as.character(subject), ...)
}

setMethod(srdistance, c("DNAStringSet", "DNAStringSet"),
          .srdistance_DNAStringSet_DNAStringSet)

## tables

.stringset_tables <- function(x, n=50, ...) {
    if (length(x) == 0) {
        return(list(top=integer(0),
                    distribution=data.frame(
                      nOccurrences=integer(0),
                      nReads=integer(0))))
    }
    ## FIXME: two sorts
    srt <- srsort(x)
    r <- srrank(x)
    t <- tabulate(r)
    o <- order(t, decreasing=TRUE)
    ## n most common sequences
    n <- min(n, sum(t!=0))              # remove duplicates
    top <- head(t[o], n)
    names(top) <- as.character(head(srt[o], n))
    ## overall frequency -- equivalent of table(table(sread))
    tt <- tabulate(t)
    nOccurrences <- seq_along(tt)[tt!=0]
    nReads <- tt[tt!=0]
    ## results
    list(top=top,
         distribution=data.frame(
           nOccurrences=nOccurrences,
           nReads=nReads, row.names=NULL))
}

setMethod(tables, "XStringSet", .stringset_tables)

## trimTails

setMethod(trimTails, "BStringSet",
    function(object, k, a, successive=FALSE, ..., alphabet,
             ranges=FALSE)
{
    if (missing(alphabet))
        alphabet <- sapply(as.raw(0:127), rawToChar)
    tryCatch({
        k <- as.integer(k)
        if (1L != length(k) || k < 0L)
            stop("'k' must be integer(1) >= 0L")
        a <- as.character(a)
        if (1L != length(a) || 1L != nchar("A"))
            stop("'", a, "' must satsify 'nchar(a) == 1L'")
        if (!a %in% alphabet)
            stop("'", a, "' must be a character with encoding < 128")
        successive <- as.logical(successive)
        if (1L != length(successive) || is.na(successive))
            stop("'successive' must be logical(1), not NA")
    }, error=function(err) {
        .throw(SRError("UserArgumentMismatch", conditionMessage(err)))
    })
    tryCatch({
        a_map <- rev(cumsum(rev(alphabet==a))) # '1' if < a
        names(a_map) <- alphabet
        ends <- .Call(.trimTails, object, k, a_map, successive)
    }, error=function(err) {
        .throw(SRError("InternalError", conditionMessage(err)))
    })
    if (ranges) IRanges(1, ends)
    else narrow(object, 1L, ends)[0L != ends]
})

setMethod(trimTails, "XStringQuality",
    function(object, k, a, successive=FALSE, ..., ranges=FALSE)
{
    rng <- callGeneric(as(object, "BStringSet"), k, a, successive, ...,
                       ranges=TRUE)
    if (ranges) rng
    else narrow(object, 1L, end(rng))[0L != width(rng)]
})
