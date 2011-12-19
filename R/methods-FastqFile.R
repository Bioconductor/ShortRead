FastqFile <-
    function(con, ...)
{
    .ShortReadFile(.FastqFile_g, con, ...)
}

setMethod(readFastq, "FastqFile",
    function(dirPath, pattern=character(), ...)
{
    if (length(pattern) != 0)
        .throw(SRWarn("UserArgumentMismatch",
                      "'pattern' ignored for '%s'",
                      "readFastq,FastqFile-method"))
    callGeneric(path(dirPath), ...)
})

setMethod(writeFastq, c("ShortReadQ", "FastqFile"),
    function(object, file, mode="w", full=FALSE, ...)
{
    if (missing(mode))
        tryCatch(mode <- summary(file$con)$mode,
                 error=function(...) NULL)
    callGeneric(object, path(file), mode=mode, full=full, ...)
})

FastqFileList <-
    function(..., class="FastqFile")
{
    Rsamtools:::.RsamtoolsFileList(..., class=class)
}
