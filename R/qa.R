.qa_character <-
    function(dirPath, pattern=character(0),
             type=c("SolexaExport", "SolexaRealign", "Bowtie",
               "MAQMap", "MAQMapShort", "fastq"),
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
           SolexaRealign=.qa_SolexaRealign(dirPath, pattern,
             type="SolexaRealign", ...),
           Bowtie=.qa_Bowtie(dirPath, pattern, type="Bowtie", ...),
           MAQMap=.qa_MAQMap(dirPath, pattern, type=type, ...),
           MAQMapShort=.qa_MAQMap(dirPath, pattern, type=type, ...),
           fastq=.qa_fastq(dirPath, pattern, type="fastq", ...))
}

setMethod(qa, "character", .qa_character)
