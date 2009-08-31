setMethod(rbind, ".QA", function(..., deparse.level=NA) {
    lst <- list(...)
    if (length(unique(sapply(lst, class))) != 1)
        .throw(SRError("UserArgumentMismatch",
                       "rbind,.QA-method '...' arguments must all be the same class"))

    f <- function(nm, lst) {
        elts <- lapply(lst, "[[", nm)
        if (class(elts[[1]]) == "list") {
            nms <- names(elts[[1]])
            l <- lapply(nms, f, elts)
            names(l) <- nms
            l
        } else {
            do.call(rbind, elts)
        }
    }
    nms <- names(.srlist(lst[[1]]))
    l <- sapply(nms, f, lst, simplify=FALSE)
    names(l) <- nms
    new(class(lst[[1]]), .srlist=l)
})

setMethod(show, ".QA", function(object) {
    callNextMethod()
    .dims <- function(elt) {
        switch(class(elt),
               matrix=,
               data.frame=paste(dim(elt), collapse=" "),
               length(elt))
    }
    .names <- function(lst, depth=0) {
        nms <- names(lst)
        for (i in seq_along(lst)) {
            fmt <- paste("%", depth*2, "s%s: %s(%s)\n", sep="")
            cat(sprintf(fmt, "", nms[i], class(lst[[i]]), .dims(lst[[i]])))
            if (is.list(lst[[i]]) && !is.data.frame(lst[[i]]))
                .names(lst[[i]], depth=depth+1)
        }
    }
    cat("QA elements (access with qa[[\"elt\"]]):\n")
    .names(.srlist(object), depth=1)
})
