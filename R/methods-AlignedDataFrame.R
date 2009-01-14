AlignedDataFrame <- function(data, metadata, nrow=nrow(data)) {
    if (missing(data)) {
        data <- data.frame(rep(0L, nrow))[,FALSE]
        metadata <- data.frame(labelDescription=character(0))
    }
    new("AlignedDataFrame", data=data, varMetadata=metadata)
}

setMethod(append, c("AlignedDataFrame", "AlignedDataFrame", "missing"),
    function(x, values, after=length(x))
{
    if (!identical(varMetadata(x), varMetadata(values))) {
        .throw(SRError("IncompatibleTypes",
                       "'%s' and '%s' have different '%s'",
                       "x", "values", "varMetadata"))
    }
    new(class(x), data=rbind(pData(x), pData(values)),
        varMetadata=varMetadata(x))
})
