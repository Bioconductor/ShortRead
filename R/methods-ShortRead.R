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
    
## subset

setMethod("[", c("ShortRead", "missing", "missing"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("ShortRead", "missing", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

setMethod("[", c("ShortRead", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA) .subset_err())

.ShortRead_subset <- function(x, i, j, ..., drop=TRUE) {
    if (nargs() != 2) .subset_err()
    initialize(x, sread=sread(x)[i], id=id(x)[i])
}

setMethod("[", c(x="ShortRead", i="ANY", j="missing"),
          .ShortRead_subset)

## manip

.abc_ShortRead <- function(stringSet, alphabet, ...) {
    if (!missing(alphabet))
        .throw(SRWarn("UserArgumentMismatch", "'alphabet' ignored"))
    alphabetByCycle(sread(stringSet), ...)
}

setMethod("alphabetByCycle", "ShortRead", .abc_ShortRead)

setMethod("clean", "ShortRead", .forward_obj)

setMethod("srorder", "ShortRead", .forward_x)

setMethod("srrank", "ShortRead", .forward_x)

setMethod("srsort", "ShortRead", .forward_x)

setMethod("srduplicated", "ShortRead", .forward_x)

setMethod("tables", "ShortRead", .forward_x)

## show

setMethod("show", "ShortRead", function(object) {
    callNextMethod()
    wd <- width(object)
    if (length(wd)>2) wd <- paste(range(wd), collapse="..")
    cat("length:", length(object), "reads; width:", wd, "cycles\n")
})

setMethod("detail", "ShortRead", function(object, ...) {
    cat("class: ", class(object), "\n")
    cat("\nsread:\n")
    show(sread(object))
    cat("\nid:\n")
    show(id(object))
})
