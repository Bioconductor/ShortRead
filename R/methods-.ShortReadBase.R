setMethod(append,
          c(".ShortReadBase", ".ShortReadBase", "ANY"),
          function(x, values, after=length(x))
{
    .throw(SRError("UserArgumentMismatch",
                   "'%s' methods not defined for classes '%s', '%s'",
                   "append", class(x), class(values)))
})

setMethod(show,
          signature=signature(object=".ShortReadBase"),
          function(object) {
              cat("class: ", class(object), "\n", sep="")
          })

setMethod(detail,
          signature=signature(object=".ShortReadBase"),
          function(object, ...) {
              cat("class: ", class(object), "\n", sep="")
          })
