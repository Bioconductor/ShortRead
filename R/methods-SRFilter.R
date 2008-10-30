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
    if (is(x, "AlignedRead"))
        alphabetFrequency(sread(x), ...)
    else
        alphabetFrequency(x, ...)
}

chromosomeFilter <-
    function(regex=character(0), .name="ChromosomeFilter")
{
    .check_type_and_length(regex, "character", 0:1)
    srFilter(function(x) {
        .idx <- logical(length(x))
        .idx[grep(regex, chromosome(x))] <- TRUE
        .idx
    }, name=.name)
}

positionFilter <- function(min=-Inf, max=Inf,
                           .name="PositionFilter")
{
    .check_type_and_length(min, "numeric", 1)
    .check_type_and_length(max, "numeric", 1)
    srFilter(function(x) {
        position(x) >= min & position(x) <= max
    }, name=.name)
}

uniqueFilter <- function(withSread=TRUE,
                         .name="UniqueFilter")
{
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
        apply(alf[,nuc], 1, max) <= threshold
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
