AlignedDataFrame <- function(data, metadata, nrow=nrow(data)) {
    if (missing(data)) {
        data <- data.frame(rep(0L, nrow))[,FALSE]
        metadata <- data.frame(labelDescription=character(0))
    }
    new("AlignedDataFrame", data=data, varMetadata=metadata)
}
