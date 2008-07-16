.qa_character <- function(dirPath, pattern=character(0),
                          type=c("SolexaExport"), ...) {
    tryCatch(type <- match.arg(type),
             error=function(err) {
                 .throw(SRError("UserArgumentMismatch",
                                conditionMessage(err)))
             })
    switch(type,
           SolexaExport=.qa_SolexaExport(dirPath, pattern,
             type="SolexaExport", ...))
}

setMethod("qa", "character", .qa_character)
