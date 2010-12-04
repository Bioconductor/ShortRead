setMethod(.srValidity, "SRFilter", function(object) {
    msg <- NULL
    fmls <- formals(object)
    if (length(fmls) != 1 || names(fmls)[[1]] != "x") {
        msg <- c(msg, paste("'filter' must have one argument, 'x'"))
    }
    if (is.null(msg)) TRUE else msg
})

setMethod(srFilter, "missing", function(fun, name, ...) {
    name <- mkScalar(as.character(name))
    new("SRFilter", function(x) !logical(length(x)),
        name=name, ...)
})
          
setMethod(srFilter, "function",
          function(fun, name, ...) 
{
    name <- mkScalar(as.character(name))
    new("SRFilter", fun, name=name, ...)
})

setMethod(srFilter, "SRFilter", function(fun, name, ...) {
    slot(fun, ".Data")
})

setMethod(name, "SRFilter", function(x, ...) {
    slot(x, "name")
})

.getAlphabetFrequency <- function(x, ...)
{
    if (is(x, "ShortRead"))
        alphabetFrequency(sread(x), ...)
    else
        alphabetFrequency(x, ...)
}

idFilter <-
    function(regex=character(0), fixed=FALSE, exclude=FALSE,
             .name="idFilter")
{
    .check_type_and_length(regex, "character", 0:1)
    srFilter(function(x) {
        .idx <- logical(length(x))
        .idx[grep(regex, as.character(id(x)), fixed=fixed)] <- TRUE
        if (exclude) .idx <- !.idx
        .idx
    }, name = .name)
}

chromosomeFilter <-
    function(regex=character(0), fixed=FALSE, exclude=FALSE,
             .name="ChromosomeFilter")
{
    .check_type_and_length(regex, "character", 0:1)
    srFilter(function(x) {
        .idx <- logical(length(x))
        .idx[grep(regex, chromosome(x), fixed=fixed)] <- TRUE
        if (exclude) .idx <- !.idx
        .idx
    }, name=.name)
}

positionFilter <-
    function(min=-Inf, max=Inf, .name="PositionFilter")
{
    .check_type_and_length(min, "numeric", 1)
    .check_type_and_length(max, "numeric", 1)
    srFilter(function(x) {
        !is.na(position(x)) & position(x) >= min &
            position(x) <= max
    }, name=.name)
}

uniqueFilter <-
    function(withSread=TRUE, .name="UniqueFilter")
{
    msg <-
        if (withSread) "occurrenceFilter(withSread=TRUE)"
        else "occurrenceFilter"
    .Deprecated(msg, package="ShortRead")
    .check_type_and_length(withSread, "logical", 1)
    srFilter(function(x) {
        if (withSread)
            !srduplicated(x)
        else {
            !(duplicated(position(x)) & duplicated(strand(x)) &
              duplicated(chromosome(x)))
        }
    }, name=.name)
}

## withSread
##   TRUE: sread, chromosome, position, strand
##   FALSE: chromosome, position, strand
##   NA: sread
.occurrenceName <-
    function(min, max, withSread, duplicates)
{
    if (!is.character(duplicates))
    {
        duplicates <- deparse(substitute(duplicates, env=parent.frame()))
        if (length(duplicates) > 1)
            duplicates <- "custom"
    }
    sprintf("%s\n  min=%d max=%d withSread='%s'\n  duplicates='%s'",
            "OccurrenceFilter", min, max, withSread, duplicates)
}

occurrenceFilter <-
    function(min=1L, max=1L, withSread=c(TRUE, FALSE, NA),
             duplicates=c("head", "tail", "sample", "none"),
             .name=.occurrenceName(min, max, withSread,
                 duplicates))
{
    .check_type_and_length(min, "numeric", 1L)
    .check_type_and_length(max, "numeric", 1L)
    if (missing(withSread))
        withSread <- withSread[1]
    .check_type_and_length(withSread, "logical", 1L)
    if (is.character(duplicates))
        duplicates <- match.arg(duplicates)
    if (max < min)
        .throw(SRError("UserArgumentMismatch",
                       "'min' must be <= 'max'"))
    srFilter(function(x) {
        rnk <- 
            if (is(x, "AlignedRead")) {
                if (is.na(withSread)) srrank(sread(x))
                else srrank(x, withSread=withSread)
            } else srrank(x)
        t <- tabulate(rnk)
        result <- rnk %in% which(t >= min & t <= max)
        if (!(is.character(duplicates) && "none" == duplicates)) {
            q <- which(rnk %in% which(t > max))
            if(length(q) != 0L) {
                x <- tapply(q, rnk[q], duplicates, max, simplify=FALSE)
                result[unlist(x, use.names=FALSE)] <- TRUE
            }
        }
        result
    }, name=.name)
}

strandFilter <- function(strandLevels=character(0),
                         .name="StrandFilter")
{
    .check_type_and_length(strandLevels, "character", NA)
    srFilter(function(x) strand(x) %in% strandLevels,
             name=.name)
}

nFilter <- function(threshold=0L, .name="CleanNFilter") 
{
    .check_type_and_length(threshold, "numeric", 1)
    srFilter(function(x) {
        .getAlphabetFrequency(x, baseOnly=TRUE)[,"other"] <= threshold
    }, name=.name)
}

polynFilter <- function(threshold=0L,
                        nuc=c("A", "C", "T", "G", "other"),
                        .name="PolyNFilter")
{
    .check_type_and_length(threshold, "numeric", 1)
    .check_type_and_length(nuc, "character", NA)
    ok <- eval(formals()[["nuc"]])
    if (!all(nuc %in% ok))
        .arg_mismatch_value_err("nuc",
                                paste(nuc, collapse=", "),
                                ok)
    srFilter(function(x) {
        alf <- .getAlphabetFrequency(x, baseOnly=TRUE)
        apply(alf[,nuc,drop=FALSE], 1, max) <= threshold
    }, name=.name)
}

srdistanceFilter <- function(subject=character(0), threshold=0L,
                             .name="SRDistanceFilter")
{
    .check_type_and_length(subject, "character", NA)
    .check_type_and_length(threshold, "numeric", 1)
    srFilter(function(x) {
        .idx <- !logical(length(x))
        dist <- srdistance(x, subject)
        for (i in seq_along(dist))
            .idx <- .idx & dist[[i]] >= threshold
        .idx
    }, name=.name)
}

dustyFilter <-
    function(threshold=Inf, batchSize=NA, .name="DustyFilter")
{
    .check_type_and_length(threshold, "numeric", 1)
    srFilter(function(x) dustyScore(x, batchSize) <= threshold,
             name=.name)
}

alignQualityFilter <- function(threshold=0L,
                               .name="AlignQualityFilter")
{
    .check_type_and_length(threshold, "numeric", 1)
    srFilter(function(x) quality(alignQuality(x)) >= threshold,
             name=.name)
}

alignDataFilter <- function(expr=expression(),
                            .name="AlignDataFilter")
{
    .check_type_and_length(expr, "expression", NA)
    srFilter(function(x) eval(expr, pData(alignData(x))),
             name=.name)
}

compose <- function(filt, ..., .name) {
    lst <- if (missing(filt)) list(...) else list(filt, ...)
    for (`filt, ...` in lst)
        .check_type_and_length(`filt, ...`, "SRFilter", NA)
    if (missing(.name))
        .name <- paste(sapply(lst, name), collapse=" o ")
    srFilter(function(x) {
        .idx <- !logical(length(x))
        for (elt in rev(lst))
            .idx <- .idx & elt(x)
        .idx
    }, name =.name)
}

setMethod(show, "SRFilter", function(object) {
    cat("class:", class(object), "\n")
    cat("name:", name(object), "\n")
    cat("use srFilter(object) to see filter\n")
})
