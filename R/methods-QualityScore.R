## interface

## constructors, [, [[, length, width, show, detail

## QualityScore

.QualityScore_subset <- function(x, i, j, ..., drop=TRUE) {
    if (nargs() != 2) .subset_err()
    initialize(x, quality=quality(x)[i])
}

setMethod("[", c("QualityScore", "ANY", "missing"),
          .QualityScore_subset)

.QualityScore_subset2 <- function(x, i, j, ...) {
    if (nargs() != 2) .subset_err()
    quality(x)[[i]]
}

setMethod("[[", c("QualityScore", "ANY", "missing"),
          .QualityScore_subset2)

setMethod("length", "QualityScore", function(x) length(quality(x)))

setMethod("width", "QualityScore", function(x)
          .undefined_method_err(class(x), "width"))

setMethod("detail", "QualityScore", function(object) {
    callNextMethod()
    cat("quality:\n")
    print(quality(object))
})

## NumericQuality

setMethod("width", "NumericQuality",
          function(x) rep(1, length(x)))

setMethod("show", "NumericQuality", function(object) {
    callNextMethod()
    .show_some("quality", quality(object))
})

## IntegerQuality

IntegerQuality <- function(quality=integer(0)) {
    new("IntegerQuality", quality=quality)
}

## MatrixQuality

MatrixQuality <- function(quality=new("matrix")) {
    new("MatrixQuality", quality=quality)
}

.MatrixQuality_subset <- function(x, i, j, ..., drop=FALSE) {
    if (nargs() != 2) .subset_err()
    initialize(x, quality=quality(x)[i,, drop=FALSE])
}

setMethod("[", c("MatrixQuality", "ANY", "missing"),
          .MatrixQuality_subset)

.MatrixQuality_subset2 <- function(x, i, j, ...) {
    if (nargs() != 2) .subset_err()
    quality(x)[i,]
}

setMethod("[[", c("MatrixQuality", "ANY", "missing"),
          .MatrixQuality_subset2)

setMethod("length", "MatrixQuality",
          function(x) nrow(x))

setMethod("width", "MatrixQuality",
          function(x) rep(ncol(x), nrow(x)))

## FastqQuality, SFastqQuality

FastqQuality <- function(quality=BStringSet(character(0))) {
    new("FastqQuality", quality=quality)
}

SFastqQuality <- function(quality=BStringSet(character(0))) {
    new("SFastqQuality", quality=quality)
}

setAs("FastqQuality", "matrix", function(from) {
    if (!length(unique(width(from)))==1)
        .throw(SRError("UserArgumentMismatch",
                       "matrix requires identical quality score widths"))
    .Call(.alphabet_as_int, quality(from), 0:255-32L)
})

setAs("SFastqQuality", "matrix", function(from) {
    if (!length(unique(width(from)))==1)
        .throw(SRError("UserArgumentMismatch",
                       "matrix requires identical quality score widths"))
    .Call(.alphabet_as_int, quality(from), 0:255-64L)
})

setMethod("width", "FastqQuality",
          function(x) width(quality(x)))

setMethod("show", "FastqQuality", function(object) {
    callNextMethod()
    cat("quality:\n")
    show(quality(object))
})

.FastqQuality_af<- function(x, baseOnly=FALSE, freq=FALSE, ...) {
    callGeneric(quality(x), freq=freq, ...)
}

setMethod("alphabetFrequency", "FastqQuality", .FastqQuality_af)

.FastqQuality_abc<- function(stringSet, alphabet, ...) {
   if (missing(alphabet))
     .abc_BStringSet(quality(stringSet),
                     alphabet=sapply(as.raw(33:132), rawToChar), ...)
   else
     .abc_BStringSet(quality(stringSet), alphabet=alphabet, ...)
}

setMethod("alphabetByCycle", "FastqQuality", .FastqQuality_abc)

.SFastqQuality_ascore<- function(object, score=0:255-64, ...) {
    .Call(.alphabet_score, quality(object), score)
}

setMethod("alphabetScore", "SFastqQuality", .SFastqQuality_ascore)

setMethod("srrank", "FastqQuality", .forward_xq)

setMethod("srorder", "FastqQuality", .forward_xq)

setMethod("srsort", "FastqQuality", .forward_xq)

setMethod("srduplicated", "FastqQuality", .forward_xq)

.FastqQuality_srduplicated<- function(x, incomparables=FALSE, ...) {
    callGeneric(x=quality(x), ...)
}

setMethod("srduplicated", "FastqQuality", .FastqQuality_srduplicated)
