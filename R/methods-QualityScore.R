## interface

## constructors, [, [[, length, width, append, show, detail

## QualityScore

.QualityScore_subset <- function(x, i, j, ..., drop=TRUE) {
    if (!missing(...)) .subset_err()
    initialize(x, quality=quality(x)[i])
}

setMethod("[", c("QualityScore", "ANY", "missing"),
          .QualityScore_subset)

.QualityScore_subset2 <- function(x, i, j, ...) {
    if (!missing(...)) .subset_err()
    quality(x)[[i]]
}

setMethod("[[", c("QualityScore", "ANY", "missing"),
          .QualityScore_subset2)

setMethod(length, "QualityScore", function(x) length(quality(x)))

setMethod(width, "QualityScore", function(x)
          .undefined_method_err(class(x), "width"))

setMethod(append, c("QualityScore", "QualityScore"),
    function(x, values, after=length(x))
{
    initialize(x, quality=append(quality(x), quality(values)))
})

setMethod(detail, "QualityScore", function(x) {
    callNextMethod()
    cat("quality:\n")
    print(quality(x))
})

## NumericQuality

setMethod(width, "NumericQuality",
          function(x) rep(1, length(x)))

setMethod(show, "NumericQuality", function(object) {
    callNextMethod()
    .show_some("quality", quality(object))
})

## IntegerQuality

IntegerQuality <- function(quality=integer(0)) {
    new("IntegerQuality", quality=quality)
}

## Import integer qualities from 454 .qual files
.readFASTA <-
    ## from Biostrings; legacy code but handles numeric values
    function(file, checkComments=TRUE, strip.descs=TRUE)
{
    if (missing(strip.descs))
        warning("use 'strip.descs=FALSE' for compatibility with old version\n",
                "  of readFASTA(), or 'strip.descs=TRUE' to remove the \">\"\n",
                "  at the beginning of the description lines and to get\n",
                "  rid of this warning (see '?readFASTA' for more details)")
    if (is.character(file)) {
        file <- file(file, "r")
        on.exit(close(file))
    } else {
        if (!inherits(file, "connection"))
            stop("'file' must be a character string or connection")
        if (!isOpen(file)) {
            open(file, "r")
            on.exit(close(file))
        }
    }

    s1 <- scan(file=file, what="", sep="\n",
               quote="", allowEscapes=FALSE, quiet=TRUE)
    if (checkComments) {
        ##comments are supposedly lines beginning with semi-colons
        comments <- grep("^;", s1)
        if (length(comments) > 0)
            s1 <- s1[-comments]
    }
    descriptions <- which(substr(s1, 1L, 1L) == ">")
    numF <- length(descriptions)
    if (numF == 0)
        stop("no FASTA sequences found")
    dp <- descriptions + 1L
    dm <- descriptions - 1L
    end <- c(dm[-1], length(s1))
    lapply(seq_len(numF),
           function(i)
           {
               desc <- s1[descriptions[i]]
               if (strip.descs)
                   desc <- substr(desc, 2L, nchar(desc))
               if (end[i] >= dp[i]) {
                   seq <- paste(s1[dp[i]:end[i]], collapse="")
               } else {
                   warning("record \"", desc, "\" contains no sequence")
                   seq <- ""
               }
               list(desc=desc, seq=seq)
           }
    )
}

.readQual <- function(file, reads = NULL) {
  if (!is.null(reads)) { ## a lot faster if the reads are known
    nums <- scan(file, integer(0), n = sum(width(reads)),
                 comment.char = ">")
    inds <- seq_len(length(reads))
    scores <- split(nums, factor(rep(inds, width(reads)), inds))
  } else {
      qual <- .readFASTA(file, strip.descs=TRUE)
      scores <- lapply(strsplit(subListExtract(qual, "seq", TRUE), " +"),
                       as.integer)
      names(scores) <- subListExtract(qual, "desc", TRUE)
  }
  scores
}

setMethod(readQual, "character",
          function(dirPath, reads = NULL, pattern=character(),
                   sample = 1, ...) 
{
    src <- .file_names(dirPath, pattern)[sample]
    scores <- do.call(c, lapply(src, .readQual, reads))
    FastqQuality(sapply(scores, function(elt) rawToChar(as.raw(elt+33))))
})

## MatrixQuality

MatrixQuality <- function(quality=new("matrix")) {
    new("MatrixQuality", quality=quality)
}

.MatrixQuality_subset <- function(x, i, j, ..., drop=FALSE) {
    if (!missing(...)) .subset_err()
    initialize(x, quality=quality(x)[i,, drop=FALSE])
}

setMethod("[", c("MatrixQuality", "ANY", "missing"),
          .MatrixQuality_subset)

.MatrixQuality_subset2 <- function(x, i, j, ...) {
    if (!missing(...)) .subset_err()
    quality(x)[i,]
}

setMethod("[[", c("MatrixQuality", "ANY", "missing"),
          .MatrixQuality_subset2)

setMethod(dim, "MatrixQuality",
          function(x) dim(quality(x)))

setMethod(length, "MatrixQuality",
          function(x) nrow(quality(x)))

setMethod(width, "MatrixQuality",
          function(x) rep(ncol(x), nrow(x)))

## FIXME: implement this, when starts are un-equal
setMethod(narrow, "MatrixQuality",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
{
    sew <- solveUserSEW(width(x), start = start, end = end, 
                        width = width)
    if (length(unique(width(sew))) != 1)
        .throw(SRError("UserArgumentMismatch", "%s of %s must be 1",
                       "'length(unique(width()))'",
                       "solved SEW"))
    if (length(unique(start(sew))) == 1) {
        idx <- unique(start(sew)) + seq_len(unique(width(sew))) - 1
        initialize(x, quality=quality(x)[,idx])
    } else {
        .throw(SRError("UserArgumentMismatch",
                       "%s requires unequal 'start' positions",
                       "'narrow,MatrixQuality-method'"))
    }
})

setMethod(append, c("MatrixQuality", "MatrixQuality"),
    function(x, values, after=length(x))
{
    initialize(x, quality=rbind(quality(x), quality(values)))
})

## FastqQuality, SFastqQuality

.FastqQuality_missing <- function(quality, ...) {
    callGeneric(BStringSet(character(0)))
}

.FastqQuality_character <- function(quality, ...) {
    callGeneric(BStringSet(quality), ...)
}

setMethod(FastqQuality, "missing", .FastqQuality_missing)

setMethod(FastqQuality, "character", .FastqQuality_character)

setMethod(FastqQuality, "BStringSet", function(quality, ...) {
    new("FastqQuality", quality=quality)
})

setMethod(SFastqQuality, "missing", .FastqQuality_missing)

setMethod(SFastqQuality, "character", .FastqQuality_character)

setMethod(SFastqQuality, "BStringSet", function(quality, ...) {
    new("SFastqQuality", quality=quality)
})

setAs("FastqQuality", "numeric", function(from) {
    v <- as.vector(t(as(from, "matrix")))
    v[!is.na(v)]
})

setAs("FastqQuality", "matrix", function(from) {
    .Call(.alphabet_as_int, quality(from), 0:255-33L)
})

setAs("FastqQuality", "PhredQuality", function(from)
{
    as(quality(from), "PhredQuality")
})

setAs("SFastqQuality", "matrix", function(from) {
    .Call(.alphabet_as_int, quality(from), 0:255-64L)
})

setAs("SFastqQuality", "SolexaQuality", function(from)
{
    as(quality(from), "SolexaQuality")
})

setMethod(width, "FastqQuality",
    function(x) width(quality(x)))

setMethod(reverse, "FastqQuality", function(x, ...) {
    do.call(class(x), list(reverse(quality(x))))
})

setMethod(narrow, "FastqQuality",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
{
    initialize(x, quality=narrow(quality(x), start, end, width,
                    use.names))
})

setMethod(alphabet, "FastqQuality",
          function(x, ...) rawToChar(as.raw(32:125), TRUE))

setMethod(encoding, "FastqQuality",
    function(x)
{
    alf <- alphabet(x)
    x <- setNames(seq(-1, length.out=length(alf)), alf)
    x[x >= 0 & x <= 41]
})

setMethod(encoding, "SFastqQuality",
    function(x)
{
    alf <- alphabet(x)
    x <- setNames(seq(-32, length.out=length(alf)), alf)
    x[x >= -5 & x <= 41]
})

setMethod(show, "FastqQuality", function(object) {
    callNextMethod()
    cat("quality:\n")
    show(quality(object))
})

.FastqQuality_af <-
    function(x, as.prob=FALSE, ...) 
{
    res <- callGeneric(quality(x), as.prob=as.prob, ...)
    if (is(res, "matrix")) {
        res <- res[,1+32:125, drop=FALSE]
        colnames(res) <- alphabet(x)
    } else {
        res <- res[1+32:125]
        names(res) <- alphabet(x)
    }
    res
}

setMethod(alphabetFrequency, "FastqQuality", .FastqQuality_af)

.FastqQuality_abc <- function(stringSet, alphabet, ...) {
    if (missing(alphabet))
        alphabet <- Biostrings::alphabet(stringSet)
   .abc_BStringSet(quality(stringSet), alphabet=alphabet, ...)
}

setMethod(alphabetByCycle, "FastqQuality", .FastqQuality_abc)

.SFastqQuality_ascore <- function(object, score=0:255-64L, ...) {
    .Call(.alphabet_score, quality(object), as.numeric(score))
}

setMethod(alphabetScore, "SFastqQuality", .SFastqQuality_ascore)

.FastqQuality_ascore <- function(object, score=0:255-33L, ...) {
    .Call(.alphabet_score, quality(object), as.numeric(score))
}

setMethod(alphabetScore, "FastqQuality", .FastqQuality_ascore)

setMethod(alphabetScore, "PhredQuality", function(object, ...) {
    .Call(.alphabet_score, object, as.numeric(0:255 - 33L))
})

setMethod(trimTailw, "FastqQuality",
    function(object, k, a, halfwidth, ..., ranges=FALSE)
{
    rng <- callGeneric(quality(object), k, a, halfwidth, ...,
                       alphabet=alphabet(object), ranges=TRUE)
    if (ranges) rng
    else narrow(object, 1L, end(rng))[0L != width(rng)]
})

setMethod(trimTails, "FastqQuality",
    function(object, k, a, successive=FALSE, ..., ranges=FALSE)
{
    rng <- callGeneric(quality(object), k, a, successive, ...,
                       alphabet=alphabet(object), ranges=TRUE)
    if (ranges) rng
    else narrow(object, 1L, end(rng))[0L != width(rng)]
})

setMethod(trimEnds, "FastqQuality",
    function(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
             ..., ranges=FALSE)
{
    rng <- callGeneric(quality(object), a, left, right, relation,
                       ..., alphabet=alphabet(object), ranges=TRUE)
    if (ranges) rng
    else narrow(object, 1L, end(rng))[0L != width(rng)]
})


setMethod(srrank, "FastqQuality", .forward_xq)

setMethod(srorder, "FastqQuality", .forward_xq)

setMethod(srsort, "FastqQuality", .forward_xq)

setMethod(srduplicated, "FastqQuality", .forward_xq)

.FastqQuality_srduplicated<- function(x, incomparables=FALSE, ...) {
    callGeneric(x=quality(x), ...)
}

setMethod(srduplicated, "FastqQuality", .FastqQuality_srduplicated)
