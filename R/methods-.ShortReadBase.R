setMethod("show",
          signature=signature(object=".ShortReadBase"),
          function(object) {
              cat("class: ", class(object), "\n", sep="")
          })

setMethod("detail",
          signature=signature(object=".ShortReadBase"),
          function(object, ...) {
              cat("class: ", class(object), "\n", sep="")
          })
