RochePath <- function(experimentPath=NA_character_,
                      readPath=experimentPath,
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
                                 run = 1, ...,
                                 nrec=-1L, skip=0L) {
  dirPath <- .file_names(readPath(dirPath)[run], pattern)[sample]
  if (any(is.na(dirPath)))
    .throw(SRError("Input/Output", "'%s' is 'NA' in '%s'",
                   "readPath", "dirPath"))
  callGeneric(dirPath, ..., nrec=nrec, skip=skip)
}

setMethod(readFasta, "RochePath", .readFasta_RochePath)

.readQual_RochePath <- function(dirPath,
                                reads = NULL,
                                pattern = "\\.qual$",
                                sample = 1,
                                run = 1, ...) {
  dirPath <- .file_names(qualPath(dirPath)[run], pattern)[sample]
  if (any(is.na(dirPath)))
    .throw(SRError("Input/Output", "'%s' is 'NA' in '%s'",
                   "qualPath", "dirPath"))
  callGeneric(dirPath, ..., reads = reads)
}

setMethod(readQual, "RochePath", .readQual_RochePath)

.readFastaQual_RochePath <- function(dirPath, fastaPattern = "\\.fna$",
                                     qualPattern = "\\.qual$", sample = 1,
                                     run = 1)
{
  reads <- readFasta(dirPath, fastaPattern, sample, run)
  quals <- readQual(dirPath, reads, qualPattern, sample, run)
  ## combine the two
  new("ShortReadQ", reads, quality=quals)
}

setMethod(read454, "RochePath",
          function(dirPath, ...) readFastaQual(dirPath, ...))

setMethod(readFastaQual, "RochePath", .readFastaQual_RochePath)

setMethod(readBaseQuality, "RochePath",
          function(dirPath, ...) .readFastaQual_RochePath(dirPath, ...))

.readFastaQual_character <- function(dirPath, fastaPattern = "\\.fna$",
                                     qualPattern = "\\.qual$", sample = 1,
                                     run = 1)
{
  callGeneric(RochePath(dirPath), fastaPattern, qualPattern, sample, run)
}

setMethod(readFastaQual, "character", .readFastaQual_character)

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
