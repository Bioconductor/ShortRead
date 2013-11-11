.qa_character <-
    function(dirPath, pattern=character(0),
             type=c("fastq", "BAM", "SolexaExport", "SolexaRealign",
               "Bowtie", "MAQMap", "MAQMapShort"),
             ...)
{
    tryCatch(type <- match.arg(type),
             error=function(err) {
                 .throw(SRError("UserArgumentMismatch",
                                conditionMessage(err)))
             })
    switch(type,
           SolexaExport=.qa_SolexaExport(dirPath, pattern, 
               type=type, ...),
           SolexaRealign=.qa_SolexaRealign(dirPath, pattern,
               type=type, ...),
           Bowtie=.qa_Bowtie(dirPath, pattern, type=type, ...),
           MAQMap=.qa_MAQMap(dirPath, pattern, type=type, ...),
           MAQMapShort=.qa_MAQMap(dirPath, pattern, type=type, ...),
           fastq=.qa_fastq(dirPath, pattern, type=type, ...),
           BAM=.qa_BAM(dirPath, pattern, type=type, ...))
}

setMethod(qa, "character", .qa_character)

setMethod(qa, "list", function(dirPath, ...)
{
    if (length(unique(sapply(dirPath, class))) != 1)
        .throw(SRError("UserArgumentMismatch",
                       "qa,list-method 'dirPath' arguments must all be of same class"))
    l <- mapply(qa, dirPath, names(dirPath), MoreArgs=list(...),
                SIMPLIFY=FALSE)
    do.call(rbind, l)
})
