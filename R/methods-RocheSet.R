.RocheSet_RochePath <- function(path, phenoData, ...) {
  if (missing(phenoData)) {
    samples <- sampleNames(path)
    runs <- runNames(path)
    df <- data.frame(samples, run=runs, row.names=1)
    phenoData <-
      new("AnnotatedDataFrame", data=df,
          varMetadata=data.frame(
            labelDescription=c("Names of sequencing runs")),
          dimLabels=c("sampleNames", "sampleColumns"))
  } else {
    if (!is(phenoData, "AnnotatedDataFrame")) {
      cls <- paste(class(phenoData), collapse=" ")
      .throw(SRError("UserArgumentMismatch",
                     "expected '%s' as '%s', but got '%s'",
                     "AnnotatedDataFrame",
                     "phenoData", cls))
    }
    dimLabels(phenoData) <- c("sampleNames", "sampleColumns")
  }

  new("RocheSet", ..., sourcePath=path, phenoData=phenoData)
}

setMethod("RocheSet", "RochePath", .RocheSet_RochePath)

setMethod("RocheSet", "character", function(path, ...) {
  .RocheSet_RochePath(RochePath(path), ...)
})
