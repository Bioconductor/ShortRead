.readBaseQuality_Solexa <-
    function(dirPath, seqPattern=character(0), prbPattern=character(0),
             ...)
{
    prbs <- readPrb(dirPath, pattern=prbPattern, ...)
	new("ShortReadQ",
        quality=prbs,
        sread =
        readXStringColumns(dirPath, pattern=seqPattern,
                           colClasses=c(rep(list(NULL), 4),
                           list("DNAString")))[[1]],
        id=BStringSet(as.character(seq_len(length(prbs)))))
}

.readBaseQuality_character <-
    function(dirPath, seqPattern=character(0), prbPattern=character(0),
             type=c(
               "Solexa"),
             ...)
{
    if (missing(type))
        .arg_missing_err("type", "readBaseQuality,character-method",
                       "help(\"readBaseQuality,character-method\")")
    if (!is.character(type) || length(type) != 1)
        .arg_mismatch_type_err("type", "character(1)")
    vals <- eval(formals(sys.function())$type)
    if (!type %in% vals)
        .arg_mismatch_value_err("type", type, vals)
    tryCatch({
        switch(type,
               Solexa=.readBaseQuality_Solexa(dirPath,
                 seqPattern=seqPattern, prbPattern=prbPattern,...))
        }, error=function(err) {
            if (is(err, "SRError")) stop(err)
            else {
                seqpat <- paste(seqPattern, collapse=" ")
                prbpat <- paste(prbPattern, collapse=" ")
                txt <- paste("'%s' failed to parse files",
                             "dirPath: '%s'",
                             "seqPattern: '%s'",
                             "prbPattern: '%s'",
                             "type: '%s'",
                             "error: %s", sep="\n  ")
                msg <- sprintf(txt, "readBaseQuality", dirPath,
                               seqpat, prbpat, type,
                               conditionMessage(err))
                .throw(SRError("Input/Output", msg))
            }
        })
}

setMethod("readBaseQuality", "character", .readBaseQuality_character)
