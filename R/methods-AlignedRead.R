###

.AlignedRead_validity <- function(object)
{
    msg <- NULL
    len <- length(sread(object))
    slts <- c("chromosome", "position", "strand", "alignQuality")
    olen <- c(length(chromosome(object)),
              length(position(object)),
              length(strand(object)),
              length(alignQuality(object)))
    if (!all(olen==len)) {
        bad <- olen!=len
        msg <- c(msg,
                 sprintf("length mismatch: expected %d, found:\n  %s",
                         len, paste(slts[bad], olen[bad], sep="=",
                                    collapse=", ")))
    }
    if (is.null(msg)) TRUE else msg
}

setMethod(.srValidity, "AlignedRead", .AlignedRead_validity)

.AlignedRead_QualityConstructor <-
    function(sread)
{
    if (length(sread) > 0)
        unlist(lapply(width(sread), polyn, nucleotides="!"))
    else
        character(0)
}

AlignedRead <- function(sread = DNAStringSet(character(0)),
                        id = BStringSet(character(length(sread))),
                        quality = FastqQuality(
                          .AlignedRead_QualityConstructor(sread)),
                        chromosome = factor(rep(NA, length(sread))),
                        position = rep(NA_integer_, length(sread)),
                        strand = factor(rep(NA_integer_, length(sread)),
                          levels=.STRAND_LEVELS),
                        alignQuality = NumericQuality(
                          rep(NA_real_, length(sread))),
                        alignData = AlignedDataFrame(
                          nrow=length(sread)))
{
    new("AlignedRead", sread=sread, id=id, quality=quality,
        chromosome=as.factor(chromosome), position=position,
        strand=strand, alignQuality=alignQuality, alignData=alignData)
}

.make_getter(c("alignQuality", "alignData"))

setMethod(chromosome, "AlignedRead",
          function(object, ...) slot(object, "chromosome"))

setMethod(position, "AlignedRead",
          function(object, ...) slot(object, "position"))

setMethod(strand, "AlignedRead",
          function(x) slot(x, "strand"))

## coerce

setAs("PairwiseAlignedXStringSet", "AlignedRead",
      function(from, to)
{
    pat <- pattern(from)
    quality <- character()
    if (is(pat, "QualityAlignedXStringSet"))
        quality <- quality(pat)
    new("AlignedRead", sread = unaligned(pat), id = names(pat),
        quality = FastqQuality(quality),
        position = start(Views(pat)),
        alignQuality = IntegerQuality(score(from)))
})

setAs("AlignedRead", "RangesList", function(from)
{
    chr <- chromosome(from)
    pos <- position(from)
    wd <- width(from)
    notNA <- !(is.na(chr) | is.na(pos) | is.na(wd))
    split(IRanges(start=pos[notNA], width=wd[notNA]), chr[notNA])
})

setAs("AlignedRead", "RangedData", function(from)
{
  chr <- chromosome(from)
  pos <- position(from)
  wd <- width(from)
  std <- strand(from)
  notNA <- !(is.na(chr) | is.na(pos) | is.na(wd) | is.na(std))
  GRanges(IRanges(pos[notNA], width=wd[notNA]), space = chr[notNA],
          id = id(from)[notNA], strand = std[notNA],
          pData(alignData(from))[notNA,,drop=FALSE])
})

setAs("AlignedRead", "GRanges", function(from)
{
    chr <- chromosome(from)
    pos <- position(from)
    wd <- width(from)
    std <- strand(from)
    notNA <- !(is.na(chr) | is.na(pos) | is.na(wd) | is.na(std))
    GRanges(chr[notNA],
            IRanges(pos[notNA], width=wd[notNA]), std[notNA],
            id = id(from)[notNA], pData(alignData(from))[notNA,,drop=FALSE])
})

setAs("AlignedRead", "GappedAlignments", function(from)
{
    if (length(from) == 0L)
        cigar <- character(0)
    else
        cigar <- paste(width(from), "M", sep="")
    GappedAlignments(rname=chromosome(from), pos=position(from),
                     cigar=cigar, strand=strand(from))
})

## subset

setMethod("[", c("AlignedRead", "missing", "missing"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("AlignedRead", "missing", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("AlignedRead", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

.AlignedRead_subset <- function(x, i, j, ..., drop=TRUE)
{
    if (0L != length(list(...))) .subset_err()
    initialize(x, sread=sread(x)[i], id=id(x)[i],
               quality=quality(x)[i],
               chromosome=factor(chromosome(x)[i]),
               position=position(x)[i], strand=strand(x)[i],
               alignQuality=alignQuality(x)[i],
               alignData=alignData(x)[i,])
}

setMethod("[", c("AlignedRead", "ANY", "missing"), .AlignedRead_subset)

setMethod(append, c("AlignedRead", "AlignedRead", "missing"),
    function(x, values, after=length(x))
{
    initialize(x,
               chromosome=.append.factor(chromosome(x),
                 chromosome(values)),
               position=append(position(x), position(values)),
               strand=.append.factor(strand(x), strand(values)),
               alignQuality=append(alignQuality(x),
                 alignQuality(values)),
               alignData=append(alignData(x), alignData(values)),
               quality=append(quality(x), quality(values)),
               sread=append(sread(x), sread(values)),
               id=append(id(x), id(values)))
})

## match, %in%

setMethod("%in%", c("AlignedRead", "RangesList"),
    function(x, table)
{
    ## could use as(x, "RangesList"), but the assumptions here (about
    ## the definition of notNA, and about split() preserving order)
    ## make this fragile enough as it is
    ##
    ## consider only sensible alignemnts
    chr <- chromosome(x)
    pos <- position(x)
    wd <- width(x)
    notNA <- !(is.na(chr) | is.na(pos) | is.na(wd))
    chr <- chr[notNA]
    ## find overlap
    rl <- split(IRanges(start=pos[notNA], width=wd[notNA]), chr)
    olap <- rl %in% table
    ## map to original indicies
    len <- seq_len(length(x))
    idx <- unlist(split(len[notNA], chr), use.names=FALSE)
    len %in% idx[unlist(olap)]
})


## srorder, etc; srsort picked up by generic

setMethod(srorder, "AlignedRead",
          function(x, ..., withSread=TRUE)
{
    if (withSread)
        order(chromosome(x), strand(x), position(x), srorder(sread(x)))
    else
        order(chromosome(x), strand(x), position(x))
})

setMethod(srrank, "AlignedRead",
          function(x, ..., withSread=TRUE)
{
    .check_type_and_length(withSread, "logical", 1)
    if (is.na(withSread))
        .throw(SRError("UserArgumentMismatch",
                       "'%s' must not be NA", "withSread"))
    o <- srorder(x)
    .Call(.aligned_read_rank, x, o, withSread, environment())
})

setMethod(srduplicated, "AlignedRead",
          function(x, ..., withSread=TRUE)
{
    duplicated(srrank(x, ..., withSread=withSread))
})

## coverage

setMethod(coverage, "AlignedRead",
    function(x, shift=0L, width=NULL, weight=1L, ...,
             coords=c("leftmost", "fiveprime"),
             extend=0L)
{
    ## Argument checking:

    if(all(is.na(chromosome(x)) == TRUE)) {
        .throw(SRError("UserArgumentMismatch",
                       "chromosome names are all 'NA' see %s",
                       '?"AlignedRead-class"'))
    }

    chrlvls <- levels(chromosome(x))
    if (!identical(shift, 0L)) {
        if (!is.numeric(shift)) {
            .throw(SRError("UserArgumentMismatch",
                           "if '%s' is not 0L, then it must be a vector of integer values\n  see %s",
                           "shift", '?"AlignedRead-class"'))
        }
        if (!all(chrlvls %in% names(shift))) {
            .throw(SRError("UserArgumentMismatch",
                           "'names(%s)' (or 'names(%s)') mismatch with 'levels(chromosome(x))'\n  see %s",
                           "shift", "start", '?"AlignedRead-class"'))
        }
        if (any(duplicated(names(shift)))) {
            .throw(SRError("UserArgumentMismatch",
                           "'names(%s)' (or 'names(%s)') have duplicates\n  see %s",
                           "shift", "start", '?"AlignedRead-class"'))
        }
    }
    if (!is.null(width)) {
        if (!is.numeric(width)) {
            .throw(SRError("UserArgumentMismatch",
                           "if '%s' is not NULL, then it must be a vector of integer values\n  see %s",
                           "width", '?"AlignedRead-class"'))
        }
        if (!all(chrlvls %in% names(width))) {
            .throw(SRError("UserArgumentMismatch",
                           "'names(%s)' (or 'names(%s)') mismatch with 'levels(chromosome(x))'\n  see %s",
                           "width", "end", '?"AlignedRead-class"'))
        }
        if (any(duplicated(names(width)))) {
            .throw(SRError("UserArgumentMismatch",
                           "'names(%s)' (or 'names(%s)') have duplicates\n  see %s",
                           "width", "end", '?"AlignedRead-class"'))
        }
    }
    if (!identical(weight, 1L)) {
        .throw(SRError("UserArgumentMismatch",
                       "weighting the reads is not supported yet, sorry\n  see %s",
                       '?"AlignedRead-class"'))
    }
    tryCatch(coords <- match.arg(coords),
        error=function(err) {
            vals <- formals(sys.function(sys.parent(4)))[["cvg"]]
            .throw(SRError("UserArgumentMismatch",
                           "'%s' must be one of '%s'\n  see %s", "coords",
                           paste(eval(vals), collapse="' '"),
                           '?"AlignedRead-class"'))
        })
    if (!is.integer(extend) ||
        !(length(extend) == 1 || length(extend) == length(x))) {
        .throw(SRError("UserArgumentMismatch",
                       "'%s' must be '%s'", "extend",
                       "integer(n)', n == 1 or length(x)"))
    }
    ## end of argument checking.

    if (coords == "leftmost") {
        rstart <- position(x) -
            ifelse(strand(x) == "+", 0L, extend)
        rend <- position(x) + width(x) - 1L +
            ifelse(strand(x) == "+", extend, 0L)
    } else {
        rstart <- position(x) -
            ifelse(strand(x) == "+", 0L, width(x) + extend - 1L)
        rend <- position(x) +
            ifelse(strand(x) == "+", width(x) + extend - 1L, 0L)
    }
    cvg <- RleList(lapply(structure(chrlvls, names = chrlvls),
                          function(chr, ...) {
                              idx <- chromosome(x) == chr
                              chr_rstart <- rstart[idx]
                              chr_rend <- rend[idx]
                              if (identical(shift, 0L))
                                  chr_shift <- 0L
                              else
                                  chr_shift <- shift[chr]
                              if (is.null(width))
                                  chr_width <- NULL
                              else
                                  chr_width <- width[chr]
                              coverage(IRanges(chr_rstart, chr_rend),
                                       shift=chr_shift, width=chr_width, ...)
                          },
                          ...),
                   compress = FALSE)
    metadata(cvg) <-
      list(method="coverage,AlignedRead-method", coords=coords, extend=extend)
    cvg
})

## show

setMethod(show, "AlignedRead", function(object) {
    callNextMethod()
    cat("chromosome:", selectSome(chromosome(object)), "\n")
    cat("position:", selectSome(position(object)), "\n")
    cat("strand:", selectSome(strand(object)), "\n")
    cat("alignQuality:", class(alignQuality(object)), "\n")
    cat("alignData varLabels:",
        selectSome(varLabels(alignData(object))), "\n")
})

setMethod(detail, "AlignedRead", function(x, ...) {
    callNextMethod()
    cat("\nchromosome:", selectSome(chromosome(x)), "\n")
    cat("position:", selectSome(position(x)), "\n")
    cat("strand:", selectSome(strand(x)), "\n")
    cat("alignQuality:\n")
    detail(alignQuality(x))
    cat("\nalignData:\n")
    show(alignData(x))
})

## summary

## perhaps summary statistics like ShortReadQ except broken down by chromosome,
## strand, and their combination
