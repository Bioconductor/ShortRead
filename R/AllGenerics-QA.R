setGeneric(".filter",
           function(object, useFilter, ...) standardGeneric(".filter"),
           signature="object")

setGeneric(".clone", function(object, ...) standardGeneric(".clone"))

setGeneric("QACollate", function(src, ...)
           standardGeneric("QACollate"))

setGeneric("qa2", function(object, state, ..., verbose=FALSE)
           standardGeneric("qa2"),
           signature="object")

setGeneric("flag", function(object, ..., verbose=FALSE)
           standardGeneric("flag"),
           signature="object")
