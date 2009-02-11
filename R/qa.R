.qa_character <-
    function(dirPath, pattern=character(0),
             type=c("SolexaExport", "MAQMapShort", "fastq"),
             ...)
{
    tryCatch(type <- match.arg(type),
             error=function(err) {
                 .throw(SRError("UserArgumentMismatch",
                                conditionMessage(err)))
             })
    switch(type,
           SolexaExport=.qa_SolexaExport(dirPath, pattern,
             type="SolexaExport", ...),
           MAQMapShort=.qa_MAQMapShort(dirPath, pattern,
             type="MAQMapShort", ...),
           fastq=.qa_fastq(dirPath, pattern, type="fastq", ...))
}

setMethod(qa, "character", .qa_character)
