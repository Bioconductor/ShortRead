RochePath <- function(basePath,
                      readPath=.srPath(basePath, "^run"),
                      qualPath=readPath,
                      ..., verbose=FALSE) {
    if (verbose) {
      .checkPath(basePath)
      .checkPath(readPath)
      .checkPath(qualPath)
    }
    new("RochePath", ..., basePath=basePath,
        readPath=readPath, qualPath=qualPath)
}

.make_getter(c("readPath", "qualPath"))

.readFasta_RochePath <- function(dirPath,
                                 pattern = "\\.fna$",
                                 sample = 1,
                                 run = 1,
                                 ...) {
  dirPath <- readPath(dirPath)[[run]]
  if (is.na(dirPath))
    .throw(SRError("Input/Output", "'%s' is 'NA' in '%s'",
                   "readPath", "dirPath"))
  callGeneric(dirPath, ..., pattern = pattern, sample = sample)
}

setMethod("readFasta", "RochePath", .readFasta_RochePath)

.readQual_RochePath <- function(dirPath,
                                pattern = "\\.qual$",
                                reads = NULL,
                                sample = 1,
                                run = 1,
                                ...) {
  dirPath <- qualPath(dirPath)[[run]]
  if (is.na(dirPath))
    .throw(SRError("Input/Output", "'%s' is 'NA' in '%s'",
                   "qualPath", "dirPath"))
  callGeneric(dirPath, ..., reads = reads, pattern = pattern, sample = sample)
}

setMethod("readQual", "RochePath", .readQual_RochePath)

.sread_RochePath <- function(object, sample = 1, run = 1) {
  reads <- readFasta(object, sample = sample, run = run)
  quals <- readQual(object, reads = reads, sample = sample, run = run)
  ## combine the two
  new("ShortReadQ", reads, quality=quals)
}

setMethod("sread", "RochePath", .sread_RochePath)

.sampleNames_RochePath <- function(object) {
  sub("_.*", "", basename(.file_names(readPath(object), "\\.fna")))
}

setMethod("sampleNames", "RochePath", .sampleNames_RochePath)

.runNames_RochePath <- function(object) {
  basename(readPath(object))
}

setMethod("runNames", "RochePath", .runNames_RochePath)
