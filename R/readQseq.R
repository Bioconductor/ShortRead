.readQseq_ShortReadQ <- 
    function(dirPath, pattern=character(0), ...,
             filtered=FALSE, verbose=FALSE) 
{
    colClasses <- rep(list(NULL), 11)
    colClasses[9:10] <- c("DNAString", "BString")
    elts <- readXStringColumns(dirPath, pattern, colClasses, ...)
    if (filtered) {
        what <- rep(list(NULL), 11)
        what[[11]] <- integer(0)
        filt <- sapply(.file_names(dirPath, pattern), function(fl, ...) {
            scan(fl, ...)[[11]] == 1
        }, what=what, quiet=!verbose)
        elts[[1]] <- elts[[1]][filt]
        elts[[2]] <- elts[[2]][filt]
    }
    new("ShortReadQ", ..., sread=elts[[1]], 
        quality=SFastqQuality(elts[[2]]),
        id=BStringSet(rep("", length(elts[[1]]))))
}

.readQseq_DataFrame <- 
    function(dirPath, pattern=character(0), ...,
             what=list(machine=character(0), run=integer(0),
               lane=integer(0), tile=integer(0), x=integer(0),
               y=integer(0), index=integer(0), readNumber=integer(0),
               sread=DNAStringSet(character(0)),
               quality=BStringSet(character(0)),
               filter=factor(levels=c("N", "Y"))),
             filtered=FALSE,
             verbose=FALSE)
{
    if (!is.list(what) || length(what) != 11)
        .arg_mismatch_type_err("what", "list(1)")
    xWhat <- what
    xstrings <- which(sapply(what, class) %in%
                      c("DNAStringSet", "BStringSet"))
    what[xstrings] <- list(NULL)
    fls <- .file_names(dirPath, pattern)
    elts <- lapply(fls, scan, what, ..., quiet=!verbose)
    data <- do.call(mapply, c(c, elts))
    if (length(xstrings) != 0) {
        xWhat[-xstrings] <- list(NULL)
        xWhat[xstrings] <- 
            lapply(xWhat[xstrings],
                   function(elt) sub("Set$", "", class(elt)))
        data[xstrings] <- readXStringColumns(dirPath, pattern, xWhat)
    }
    xdf <- do.call(DataFrame, data)
    if (is.factor(what[[11]])) {
        xdf[[11]] <- factor(levels(what[[11]])[xdf[[11]] + 1],
                            levels=levels(what[[11]]))
        if (filtered)
            xdf <- xdf[xdf[[11]] == "Y", -11]
    }
    xdf
}

.readQseq_character <-
    function(dirPath, pattern=character(0), ...,
             as=c("ShortReadQ", "DataFrame", "XDataFrame"),
             filtered=FALSE,
             verbose=FALSE)
{
    if (missing(as)) {
        as <- "ShortReadQ"
    } else if (!is.character(as) || length(as) != 1) {
        .arg_mismatch_type_err("as", "character(1)")
    } else {
        vals <- eval(formals(ShortRead:::.readQseq_character)$as)
        if (!as %in% vals)
            .arg_mismatch_value_err("as", as, vals)
    }

    tryCatch({
        switch(as,
               ShortReadQ=.readQseq_ShortReadQ(
                 dirPath, pattern, ...,
                 filtered=filtered, verbose=verbose),
               DataFrame=.readQseq_DataFrame(
                 dirPath, pattern, ...,
                 filtered=filtered, verbose=verbose),
               XDataFrame={
                   .Defunct(msg="Use type='DataFrame' instead")
               })
    }, error=function(err) {
        if (is(err, "SRError")) stop(err)
        else {
            txt <- paste("'%s' failed to parse files",
                         "dirPath: '%s'",
                         "pattern: '%s'",
                         "as: '%s'",
                         "error: %s", sep="\n  ")
            msg <- sprintf(txt, "readQseq",
                           paste(dirPath, collapse="'\n    '"),
                           pattern, as, conditionMessage(err))
            .throw(SRError("Input/Output", msg))
        }
    })
}

setMethod(readQseq, "character", .readQseq_character)
