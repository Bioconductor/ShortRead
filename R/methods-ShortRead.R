.ShortRead_validity <- function(object) {
    msg <- NULL
    if (length(sread(object)) != length(id(object)))
        msg <- c(msg,
                 sprintf("sread() and id() length mismatch: %d, %d",
                         length(sread(object)), length(id(object))))
    if (is.null(msg)) TRUE else msg
}

setMethod(".srValidity", "ShortRead", .ShortRead_validity)

.make_getter(slotNames("ShortRead"))

setMethod("length", "ShortRead", function(x) length(sread(x)))

setMethod("width", "ShortRead", function(x) {
    if (length(sread(x)) > 0) {
        unique(width(sread(x)))
    } else {
        0
    }
})
    
.ShortRead_subset_err <- function() {
    .throw(SRError("UserSubset",
                  "'[' must be called with only subscript 'i'"))
}

setMethod("[", c("ShortRead", "missing", "missing"),
          function(x, i, j, ..., drop=NA) .ShortRead_subset_err())

setMethod("[", c("ShortRead", "missing", "ANY"),
          function(x, i, j, ..., drop=NA) .ShortRead_subset_err())

setMethod("[", c("ShortRead", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA) .ShortRead_subset_err())

.ShortRead_subset <- function(x, i, j, ..., drop=TRUE) {
    if (nargs() != 2) .ShortRead_subset_err()
    initialize(x, sread=sread(x)[i], id=id(x)[i])
}

setMethod("[", c(x="ShortRead", i="ANY", j="missing"),
          .ShortRead_subset)

setMethod("show", "ShortRead", function(object) {
    callNextMethod()
    cat("length:", length(object), "reads; width:", width(object),
    "cycles\n")
})

setMethod("detail", "ShortRead", function(object, ...) {
    cat("class: ", class(object), "\n")
    cat("\nsread:\n")
    show(sread(object))
    cat("\nid:\n")
    show(id(object))
})
