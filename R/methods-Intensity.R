## IntensityMeasure

setMethod(show, "IntensityMeasure", function(object) 
{
    callNextMethod()
    cat("  dim: ", dim(object), "\n")
})

setMethod(get("["), c("IntensityMeasure", "ANY", "ANY", "ANY"),
          function(x, i, j, ..., drop=FALSE)
{
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    initialize(x, x@.Data[i,j,...,drop=FALSE])
})

setMethod(get("[["), c("ArrayIntensity", "ANY", "ANY"),
          function(x, i, j, k, ...) 
{
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    if (missing(k)) k <- TRUE
    x@.Data[i,j,k]
})

## IntensityInfo

## Intensity

setMethod(.srValidity, "Intensity", function(object) 
{
    msg <- NULL
    if (.hasMeasurementError(object) &&
        !all(dim(intensity(object)) == dim(measurementError(object)))) {
        msg <- c(msg,
                 "'intensity' and 'measurementError' dimensions differ")
    }
    if (nrow(readInfo(object)) != nrow(intensity(object))) {
        msg <- c(msg,
                 "'intensity' and 'readInfo' read numbers differ")
    }
    if (is.null(msg)) TRUE else msg
})

measurementError <- function(object, ...)
{
    if (!.hasMeasurementError(object))
        .throw(SRError("ValueUnavailable",
                       "'%s' has no value '%s'",
                       class(object), "nse"))
    slot(object, "measurementError")
}

local({
    slts <- slotNames("Intensity")
    .make_getter(slts[slts!="measurementError"], verbose=TRUE)
})

setMethod(dim, "Intensity", function(x) 
{
    dim(intensity(x))
})

setMethod(show, "Intensity", function(object)
{
    callNextMethod()
    cat("dim:", dim(object), "\n")
    cat("readInfo:", class(readInfo(object)), "\n")
    cat("intensity:", class(intensity(object)), "\n")
    if (.hasMeasurementError(object)) {
        cat("measurementError:", class(measurementError(object)), "\n")
    } else {
        cat("measurementError: not available\n")
    }
})
