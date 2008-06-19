RochePath <- function(basePath,
                      dataPath=.srPath(basePath, "\\.fna$"),
                      qualPath=.srPath(basePath, "\\.qual$"),
                      ..., verbose=FALSE) {
    if (verbose) {
      .checkPath(basePath)
      .checkPath(dataPath)
      .checkPath(qualPath)
    }
    new("RochePath", ..., basePath=basePath,
        dataPath=dataPath, qualPath=qualPath)
}

## TODO: make accessors
