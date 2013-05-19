.filterFastq_check_fnames <- function(files, destinations)
{
    if (missing('destinations'))
        stop("'destinations' missing")
    if (length(files) != length(destinations))
        stop("'files' and 'destinations' must have same length")
    if (any(exists <- file.exists(destinations)))
        stop("'destinations' exist:\n  ",
             paste(destinations[exists], collapse="\n  "))
}

.filter1 <-
    function(filter, file, destination, ..., yieldSize)
{
    strm <- FastqStreamer(file, yieldSize)
    on.exit(close(strm))
    tot <- tot1 <- nNuc <- nNuc1 <- 0L
    while (length(fq <- yield(strm))) {
        tot <- tot + length(fq)
        nNuc <- nNuc + sum(width(fq))
        fq <- if (is(filter, "FilterRules")) {
            subsetByFilter(fq, filter)
        } else filter(fq, ...)
        tot1 <- tot1 + length(fq)
        nNuc1 <- nNuc1 + sum(width(fq))
        writeFastq(fq, destination, "a")
    }
    attr(destination, "filter") <-
        data.frame(Reads=tot, KeptReads=tot1, Nucl=nNuc, KeptNucl=nNuc1)
    destination
}

filterFastq <-
    function(files, destinations, ..., filter=FilterRules(),
             yieldSize=1000000L)
{
    if (missing(filter))
        warning("'filterFastq' invoked with missing 'filter'")
    .filterFastq_check_fnames(files, destinations)
    ## FIXME parallel over files, esp. random numbers
    x <- Map(.filter1, files, destinations,
             MoreArgs=list(filter=filter, ..., yieldSize=yieldSize))

    stats <- do.call(rbind, lapply(x, attr, "filter"))
    rownames(stats) <- make.unique(basename(files))
    attr(destinations, "filter") <- stats
    destinations
}
