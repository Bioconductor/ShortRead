ShortReadIntensityInfo <-
  function(lane=integer(0), tile=integer(length(lane)),
           x=integer(length(lane)), y=integer(length(lane)))
{
    new("AnnotatedDataFrame",
        data=data.frame(
          lane=lane, tile=tile, x=x, y=y),
        varMetadata=data.frame(
          labelDescription=c(
            "Solexa lane nubmer",
            "Solexa tile nubmer",
            "Tile x coordinate",
            "Tile y coordinate")))
}

ShortReadIntensity <-
    function(intensity=matrix(0, 0, 0), nse=matrix(0, 0, 0),
             readInfo=ShortReadIntensityInfo(
               lane=integer(nrow(intensity))),
             ...)
{
    .hasStandardErrors <- mkScalar(!missing(nse))
    new("ShortReadIntensity", intensity=intensity, nse=nse,
        readInfo=readInfo, .hasStandardErrors=.hasStandardErrors,
        ...)
}

setMethod(.srValidity, "ShortReadIntensity", function(object) 
{
    msg <- NULL
    if (.hasStandardErrors(object) &&
        !all(dim(intensity(object)) == dim(nse(object)))) {
        msg <- c(msg,
                 "'intensity' and 'nse' matrix dimensions differ")
    }
    if (nrow(readInfo(object)) != nrow(intensity(object))) {
        msg <- c(msg,
                 "'intensity' and 'readInfo' read numbers differ")
    }
    reqd <- c("lane", "tile", "x", "y")
    if (!all(reqd %in% varLabels(readInfo(object)))) {
        missing <- reqd[!reqd %in% names(readInfo(object))]
        msg <- c(msg,
                 sprintf("'readInfo' must contain columns '%s'",
                         paste(missing, collapse="' '")))
    }
    if (is.null(msg)) TRUE else msg
})

nse <- function(object, ...)
{
    if (!.hasStandardErrors(object))
        .throw(SRError("ValueUnavailable",
                       "'%s' has no value '%s'",
                       "ShortReadIntensity", "nse"))
    slot(object, "nse")
}

local({
    slts <- slotNames("ShortReadIntensity")
    .make_getter(slts[slts!="nse"], verbose=TRUE)
})

setMethod(show, "ShortReadIntensity", function(object) 
{
    callNextMethod()
    cat("readInfo:", varLabels(readInfo(object)), "\n")
    int <- intensity(object)
    cat("intensity: matrix(", nrow(int), " ", ncol(int), ")\n",
        sep="")
    if (.hasStandardErrors(object)) {
        nse <- nse(object)
        cat("nse: matrix(", nrow(nse), " ", ncol(nse), ")\n", sep="")
    } else {
        cat("nse: not available\n")
    }
})
