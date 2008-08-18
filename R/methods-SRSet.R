.SRSet_validity <- function(object) {
  msg <- NULL
  len <- length(readInd(object))
  rlen <- c(readData = nrow(readData(object)))
  if (!all(rlen==len)) {
    bad <- rlen!=len
    msg <- c(msg,
             sprintf("read length mismatch: expected %d, found:\n  %s",
                     rlen, paste(names(rlen)[bad], rlen[bad], sep="=",
                                collapse=", ")))
  }
  snames <- sampleNames(experimentPath(object))
  slen <- length(snames)
  oslen <- c(phenoData = nrow(phenoData(object)),
             readCount = length(readCount(object)))
  if (!all(oslen==slen)) {
    bad <- oslen!=slen
    msg <- c(msg,
             sprintf("sample length mismatch: expected %d, found:\n  %s",
                     slen, paste(names(oslen)[bad], oslen[bad], sep="=",
                                collapse=", ")))
  }
  osnames <- sampleNames(object)
  stest <- snames == osnames
  if (!all(stest))
    msg <- c(msg,
             sprintf("sample names mismatch:\n  %s",
                     slen, paste(snames[!stest], osnames[!stest],
                                 sep = "!=", collapse = ", ")))
  rind <- readInd(object)
  if (!all(rind > 0 & rind <= len))
    msg <- c(msg, "values in 'readInd' must be > 0 and <= number of reads")
  rcount <- readCount(object)
  if (!all(rcount >= 0))
    msg <- c(msg, "values in 'readCount' must be non-negative")
  if (sum(rcount) != len)
    msg <- c(msg,
             sprintf("'sum(readCount)', %d, must equal the number of reads, %d",
                     sum(rcount), len))
  if (is.null(msg)) TRUE else msg
}

setMethod(".srValidity", "SRSet", .SRSet_validity)

setMethod("sampleNames", "SRSet", function(object) {
  sampleNames(phenoData(object))
})
          
setMethod("show", "SRSet", function(object) {
    callNextMethod()
    cat("basePath(sourcePath(object)):\n  ",
        basePath(sourcePath(object)), "\n", sep="")
})

setMethod("detail", "SRSet", function(object, ...) {
    callNextMethod()
    cat("\n")
    detail(sourcePath(object), ...)
    cat("\nclass: AnnotatedDataFrame\n")
    pd <- phenoData(object)
    cat("pData:\n")
    print(pData(pd))
    cat("varMetadata:\n")
    print(varMetadata(pd))
})

setMethod("phenoData", "SRSet", function(object) object@phenoData)

.make_getter(c("readData", "sourcePath"))

## proposed
##setMethod("readSRQ", "SRSet", function(object) readSRQ(sourcePath(object)))
