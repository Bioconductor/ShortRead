.binReader <-
    function(con, n)
    ## read 'n' bytes from 'con', returning raw()
{
    readBin(con, raw(), n)
}

setMethod(yield, "FastqFileReader", function(x, ...) x$yield(...))
