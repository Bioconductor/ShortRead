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
setMethod("featureData", "SRSet", function(object) object@featureData)


