.readIntensities_character <-
    function(dirPath, pattern=character(0), ...,
             type=c("SolexaIntensity", "IparIntensity", "RtaIntensity"))
{
    if (missing(type)) {
        type <- "SolexaIntensity"
    } else if (!is.character(type) || length(type) != 1) {
        .arg_mismatch_type_err("type", "character(1)")
    } else {
        vals <- eval(formals(ShortRead:::.readIntensities_character)$type)
        if (!type %in% vals)
            .arg_mismatch_value_err("type", type, vals)
    }
    tryCatch({
        switch(type,
               SolexaIntensity=.readIntensities_SolexaIntensity(
                 dirPath, pattern, ...),
               IparIntensity=.readIntensities_IparIntensity(
                 dirPath, pattern, ...),
               RtaIntensity=.readIntensities_RtaIntensity(
                 dirPath, pattern, ...))
    }, error=function(err) {
        if (is(err, "SRError")) stop(err)
        else {
            pat <- paste(pattern, collapse=" ")
            txt <- paste("'%s' failed to parse files",
                         "dirPath: '%s'",
                         "pattern: '%s'",
                         "type: '%s'",
                         "error: %s", sep="\n  ")
            msg <- sprintf(txt, "readIntensities",
                           paste(dirPath, collapse="'\n    '"),
                           pat, type, conditionMessage(err))
            .throw(SRError("Input/Output", msg))
        }
    })
    
}

setMethod(readIntensities, "character", .readIntensities_character)
