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

## Import integer qualities from 454 .qual files
.readQual <- function(file, reads = NULL) {
  if (!is.null(reads)) { ## a lot faster if the reads are known
    nums <- scan(file, integer(0), n = sum(width(reads)), comment.char = ">")
    scores <- split(nums, rep(seq_len(length(reads)), width(reads)))
  } else {
    lines <- readLines(file)
    headerLines <- grep("^>", lines)
    ids <- sub("^(.*?) .*", "\\1", lines[headerLines])
    scores <- lapply(strsplit(lines[headerLines+1], " "), as.integer)
    names(scores) <- ids
  }
  scores
}
setMethod("readQual", "character", function(dirPath, reads = NULL,
                                            pattern=character(), sample = 1,
                                            ...) {
  src <- .file_names(dirPath, pattern)[sample]
  scores <- do.call("c", lapply(scores, .readQual, reads))
  IntegerQuality(scores)
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

.FastqQuality_missing <- function(quality, ...) {
    callGeneric(BStringSet(character(0)))
}

.FastqQuality_character <- function(quality, ...) {
    callGeneric(BStringSet(quality), ...)
}

setMethod("FastqQuality", "missing", .FastqQuality_missing)

setMethod("FastqQuality", "character", .FastqQuality_character)

setMethod("FastqQuality", "BStringSet", function(quality, ...) {
    new("FastqQuality", quality=quality)
})

setMethod("SFastqQuality", "missing", .FastqQuality_missing)

setMethod("SFastqQuality", "character", .FastqQuality_character)

setMethod("SFastqQuality", "BStringSet", function(quality, ...) {
    new("SFastqQuality", quality=quality)
})

setAs("FastqQuality", "numeric", function(from) {
    as.vector(as(from, "matrix"))
})

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

setMethod("alphabet", "FastqQuality",
          function(x) rawToChar(as.raw(32:125), TRUE))

setMethod("show", "FastqQuality", function(object) {
    callNextMethod()
    cat("quality:\n")
    show(quality(object))
})

.FastqQuality_af <- function(x, baseOnly=FALSE, freq=FALSE, ...) {
    res <- callGeneric(quality(x), freq=freq, ...)
    if (is(res, "matrix")) {
        res <- res[,1+32:125, drop=FALSE]
        colnames(res) <- alphabet(x)
    } else {
        res <- res[1+32:125]
        names(res) <- alphabet(x)
    }
    res
}

setMethod("alphabetFrequency", "FastqQuality", .FastqQuality_af)

.FastqQuality_abc <- function(stringSet, alphabet, ...) {
    if (missing(alphabet))
        alphabet <- Biostrings::alphabet(stringSet)
   .abc_BStringSet(quality(stringSet), alphabet=alphabet, ...)
}

setMethod("alphabetByCycle", "FastqQuality", .FastqQuality_abc)

.SFastqQuality_ascore <- function(object, score=0:255-64L, ...) {
    .Call(.alphabet_score, quality(object), as.numeric(score))
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
