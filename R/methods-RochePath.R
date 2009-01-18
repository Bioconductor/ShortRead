RochePath <- function(experimentPath=NA_character_,
                      readPath=.srPath(experimentPath, "^run"),
                      qualPath=readPath,
                      ..., verbose=FALSE) {
    if (verbose) {
      .checkPath(experimentPath)
      .checkPath(readPath)
      .checkPath(qualPath)
    }
    new("RochePath", ..., basePath=experimentPath,
        readPath=readPath, qualPath=qualPath)
}

.make_getter(c("readPath", "qualPath"))

.readFasta_RochePath <- function(dirPath,
                                 pattern = "\\.fna$",
                                 sample = 1,
                                 run, ...) {
  dirPath <- readPath(dirPath)[run]
  if (is.na(dirPath))
    .throw(SRError("Input/Output", "'%s' is 'NA' in '%s'",
                   "readPath", "dirPath"))
  callGeneric(dirPath, ..., pattern = pattern, sample = sample)
}

setMethod(readFasta, "RochePath", .readFasta_RochePath)

.readQual_RochePath <- function(dirPath,
                                pattern = "\\.qual$",
                                reads = NULL,
                                sample = 1,
                                run, ...) {
  dirPath <- qualPath(dirPath)[run]
  if (is.na(dirPath))
    .throw(SRError("Input/Output", "'%s' is 'NA' in '%s'",
                   "qualPath", "dirPath"))
  callGeneric(dirPath, ..., reads = reads, pattern = pattern, sample = sample)
}

setMethod(readQual, "RochePath", .readQual_RochePath)

.read454_RochePath <- function(dirPath, sample = 1, run = 1) {
  reads <- readFasta(dirPath, sample = sample, run = run)
  quals <- readQual(dirPath, reads = reads, sample = sample, run = run)
  ## combine the two
  new("ShortReadQ", reads, quality=quals)
}

setMethod(read454, "RochePath", .read454_RochePath)

.sampleNames_RochePath <- function(object) {
    path <- readPath(object)
    if (!is.na(path))
        sub("_.*", "", basename(.file_names(path, "\\.fna")))
    else
        callNextMethod()
}

setMethod(sampleNames, "RochePath", .sampleNames_RochePath)

.runNames_RochePath <- function(object) {
  basename(readPath(object))
}

setMethod(runNames, "RochePath", .runNames_RochePath)

setMethod(show, "RochePath", function(object) {
    callNextMethod()
    .show_additionalPathSlots(object)
})

setMethod(detail, "RochePath", function(object, ...) {
    callNextMethod()
    .detail_additionalPathSlots(object)
})
