readXStringColumns <- function(dirPath, pattern=character(0),
                               colClasses=list(NULL), sep="\t",
                               header=FALSE) {
    if (!is.list(colClasses))
        .arg_mismatch_type_err("colClasses", "list()")
    colIndex <- which(!sapply(colClasses, is.null))
    colClasses <- sub("Set$", "", colClasses[colIndex])
    okClasses <- names(slot(getClass("XString"), "subclasses"))
    if (!all(colClasses %in% okClasses)) {
        bad <- colClasses[!colClasses %in% okClasses]
        .throw(SRError("UserArgumentMismatch",
                       "'colClasses' contains invalid class%s '%s';\n  must be one of '%s'",
                       if (length(colClasses)>1) "es" else "",
                       paste(bad, collapse="' '"),
                       paste(okClasses, collapse="', '")))
    }
    files <- .file_names(dirPath, pattern)
    res <- tryCatch({
        .Call(.read_XStringSet_columns, files,
              colIndex, colClasses, sep, header)
    }, error=function(err) {
        .throw(SRError("Input/Output",
                       "while reading files '%s':\n    %s",
                       paste(basename(files), collapse=", "),
                       conditionMessage(err)))
    })
    if (header) {
        nms <- strsplit(readLines(files[[1]], 1), sep)[[1]]
        names(res) <- nms[colIndex]
    } else {
        names(res) <- names(colIndex)
    }
    res
}

