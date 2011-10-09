setMethod(.ShortReadFile, "character",
    function(g, path, open="", ...)
{
    path <- .file_names(path, character())
    g$new(con=file(path, open, encoding="ASCII"), path=path, ...)
})

setMethod(.ShortReadFile, "connection",
    function(g, path, ...)
{
    descr <- summary(path)$description
    g$new(con=path, path=descr, ...)
})

setMethod(path, "ShortReadFile", function(object, ...)
{
    object$path
})

open.ShortReadFile <-
    function(con, ...)
{
    tryCatch(open(con$con), error=function(err, ...) {
        .throw(SRError("Input/Output", "error: %s\n%s",
                       conditionMessage(err),
                       Rsamtools:::.ppath("  path", con$path)))
    })
    invisible(con)
}

close.ShortReadFile <-
    function(con, ...)
{
    tryCatch(close(con$con), error=function(err, ...) {
        .throw(SRError("Input/Output", "error: %s\n%s",
                       conditionMessage(err),
                       Rsamtools:::.ppath("  path", con$path)))
    })
    invisible(con)
}

setMethod(isOpen, "ShortReadFile", function(con, rw="") 
{
    tryCatch(isOpen(con$con), error=function(err, ...) {
        msg <- conditionMessage(err)
        if (msg != "invalid connection")
            .throw(SRWarn("Input/Output", "warning: %s\n%s",
                          conditionMessage(err),
                          Rsamtools:::.ppath("  path", con$path)))
        FALSE
    })
})
