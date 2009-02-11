.SolexaExportQA <- function(x, ...) 
{
    new("SolexaExportQA", .srlist=x, ...)
}

.report_SolexaExportQA <-
    function(x, ..., dest=paste(tempfile(), "pdf", sep="."),
             type="pdf")
{
    to <- tempfile()
    qa <- x
    save(qa, file=to)
    res <- callGeneric(to, ..., dest=dest, type=type)
    unlink(to)
    res
}

setMethod(report, "SolexaExportQA", .report_SolexaExportQA)
