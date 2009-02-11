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
