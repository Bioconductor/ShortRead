## SolexaIntensityInfo

setMethod(.srValidity, "SolexaIntensityInfo", function(object) 
{
    msg <- NULL
    reqd <- c("lane", "tile", "x", "y")
    if (slot(object, ".init")==TRUE &&
        !all(reqd %in% varLabels(object))) {
        missing <- reqd[!reqd %in% names(object)]
        msg <- c(msg,
                 sprintf("'%s' must contain columns '%s'",
                         class(object),
                         paste(missing, collapse="' '")))
    }
    if (is.null(msg)) TRUE else msg
})

## SolexaIntensity

SolexaIntensity <-
    function(intensity=array(0, c(0, 0, 0)),
             measurementError=array(0, c(0, 0, 0)),
             readInfo=SolexaIntensityInfo(
               lane=integer(nrow(intensity))),
             ...)
{
    .hasMeasurementError <- mkScalar(!missing(measurementError))
    new("SolexaIntensity",
        intensity=ArrayIntensity(intensity),
        measurementError=ArrayIntensity(measurementError),
        readInfo=readInfo,
        .hasMeasurementError=.hasMeasurementError,
        ...)
}

.readIntensities_SolexaIntensity <-
    function(dirPath, pattern=character(0), ...,
             intExtension="_int.txt",
             nseExtension="_nse.txt",
             withVariability=TRUE, verbose=FALSE)
{
    .check_type_and_length(withVariability, "logical", 1)
    .check_type_and_length(pattern, "character", NA)
    .check_type_and_length(intExtension, "character", 1)
    .check_type_and_length(nseExtension, "character", 1)
    intPattern <- paste(pattern, intExtension, sep="")
    nrec <- countLines(dirPath, intPattern)
    crec <- c(0, cumsum(nrec))
    if (withVariability) {
        nsePattern <- paste(pattern, nseExtension, sep="")
        extrec <- countLines(dirPath, nsePattern)
        if (length(nrec) != length(extrec)) {
            .throw(SRError("UserArgumentMismatch",
                           "number of files found differs between 'intensity' (%d) and 'nse' (%d)\n  dirPath: '%s'\n  pattern: '%s'\n  intExtension: '%s'\n  nseExtension: '%s'",
                           length(nrec), length(extrec),
                           dirPath, pattern, intExtension, nseExtension))
        }
        if (!all(nrec == extrec)) {
            .throw(SRError("UserArgumentMismatch",
                           "line counts differ between 'intensity' and 'nse'\n  dirPath: '%s'\n  pattern: '%s'\n  intExtension: '%s'\n  nseExtension: '%s'",
                           dirPath, pattern, intExtension, nseExtension))
        }
    }

    fls <- .file_names(dirPath, intPattern)
    gz <- gzfile(fls[[1]], "rb")
    tryCatch({
        ln <- readLines(gz, 1)
    }, finally=close(gz))
    cycles <- length(gregexpr("\t", ln, fixed=TRUE)[[1]]) - 3L
    reads <- sum(nrec)
    what <- c(rep(list(integer()), 4), rep(list(numeric()), cycles * 4L))
    int <- array(numeric(), c(reads, 4L, cycles),
                 dimnames=list(NULL, c("A", "C", "G", "T"), NULL))
    df <- data.frame(lane=integer(reads), tile=integer(reads),
                     x=integer(reads), y=integer(reads))
    for (i in seq_along(fls)) {
        tryCatch({
            gz <- gzfile(fls[[i]], "rb")
            data <- scan(gz, what, nrec[[i]],..., quiet=!verbose)
            idx <- (crec[i]+1):crec[i+1]
            int[idx,,] <- array(unlist(data[-(1:4)]),
                                c(nrec[[i]], 4L, cycles))
            df[idx,] <- data[1:4]
        }, error=function(err) {
            msg <- sprintf("parsing '%s'\n  file: %s\n  error: %s",
                           "intPattern", fls[[i]],
                           conditionMessage(err))
            .throw(SRError("Input/Output", msg))
        }, finally=close(gz))
    }
    if (withVariability) {
        fls <- .file_names(dirPath, nsePattern)
        nse <- array(numeric(), c(reads, 4L, cycles),
                     dimnames=list(NULL, c("A", "C", "G", "T"), NULL))
        what <- c(rep(list(NULL), 4), what[-(1:4)])
        for (i in seq_along(fls)) {
            tryCatch({
                gz <- gzfile(fls[[i]], "rb")
                data <- scan(gz, what, nrec[[i]], ..., quiet=!verbose)
                idx <- (crec[i]+1):crec[i+1]
                nse[idx,,] <- array(unlist(data[-(1:4)]),
                                    c(nrec[[i]], 4L, cycles))
            }, error=function(err) {
                msg <- sprintf("parsing '%s'\n  file: %s\n  error: %s",
                               "nsePattern", fls[[i]],
                               conditionMessage(err))
                .throw(SRError("Input/Output", msg))
            }, finally=close(gz))
        }
    }
    readInfo <- SolexaIntensityInfo(df[[1]], df[[2]], df[[3]], df[[4]])
    if (withVariability)
        SolexaIntensity(int, nse, readInfo)
    else
        SolexaIntensity(int, readInfo=readInfo)
}

setMethod(get("["), c("SolexaIntensity", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=TRUE)
{
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    if (missing(k)) k <- TRUE
    if (.hasMeasurementError(x))
        initialize(x, intensity=intensity(x)[i,j,k],
                   measurementError=measurementError(x)[i,j,k],
                   readInfo=readInfo(x)[i,])
    else
        initialize(x, intensity=intensity(x)[i,j,k],
                   readInfo=readInfo(x)[i,])
})
