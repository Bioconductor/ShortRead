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

setMethod("width", "FastqQuality",
          function(x) width(quality(x)))

setMethod("show", "FastqQuality", function(object) {
    callNextMethod()
    cat("quality:\n")
    show(quality(object))
})

.af_FastqQuality <- function(x, baseOnly=FALSE, freq=FALSE, ...) {
    callGeneric(quality(x), baseOnly=baseOnly, freq=freq, ...)
}

setMethod("alphabetFrequency", "FastqQuality", .af_FastqQuality)

.abc_FastqQuality <- function(stringSet, alphabet, ...) {
   if (missing(alphabet))
     .abc_BStringSet(quality(stringSet), ...)
   else
     .abc_BStringSet(quality(stringSet), alphabet=alphabet, ...)
}

setMethod("alphabetByCycle", "FastqQuality", .abc_FastqQuality)
