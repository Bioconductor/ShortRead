.readIntensities_character <-

    function(dirPath, pattern=character(0), ...,
             intExtension="_int.txt", nseExtension="_nse.txt",
             withNse=TRUE, verbose=FALSE)
{
    .check_type_and_length(withNse, "logical", 1)
    .check_type_and_length(pattern, "character", NA)
    .check_type_and_length(intExtension, "character", 1)
    .check_type_and_length(nseExtension, "character", 1)
    intPattern <- paste(pattern, intExtension, sep="")
    nrec <- countLines(dirPath, intPattern)
    crec <- c(0, cumsum(nrec))
    if (withNse) {
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
    ln <- readLines(fls[[1]], 1)
    cycles <- length(gregexpr("\t", ln, fixed=TRUE)[[1]]) - 3L
    nms <- paste(c("A", "C", "G", "T"), rep(seq_len(cycles), each=4),
                 sep="")
    ncol <- length(nms); nrow <- sum(nrec)
    what <- c(rep(list(integer()), 4), rep(list(numeric()), ncol))
    int <- matrix(numeric(), nrow=nrow, ncol=ncol,
                  dimnames=list(NULL, nms))
    df <- data.frame(lane=integer(nrow), tile=integer(nrow),
                     x=integer(nrow), y=integer(nrow))
    for (i in seq_along(fls)) {
        if (verbose)
            cat(".readIntensities_character (intPattern)\n ", fls[i],
                "\n")
        tryCatch({
            data <- scan(fls[i], what, ..., quiet=!verbose)
            idx <- (crec[i]+1):crec[i+1]
            int[idx,] <- matrix(unlist(data[-(1:4)]), ncol=ncol)
            df[idx,] <- data[1:4]
        }, error=function(err) {
            .throw(SRError("Input/Output",
                           sprintf("parsing 'intPattern'\n  file: %s\n  error: %s",
                                   fls[[i]], conditionMessage(err))))
        })
    }
    if (withNse) {
        fls <- .file_names(dirPath, nsePattern)
        nse <- matrix(numeric(), nrow=nrow, ncol=ncol,
                      dimnames=list(NULL, nms))
        for (i in seq_along(fls)) {
            if (verbose)
                cat(".readIntensities_character (nsePattern)\n ",
                    fls[i], "\n")
            tryCatch({
                what <- c(rep(list(NULL), 4), what[-(1:4)])
                data <- scan(fls[i], what, ..., quiet=!verbose)
                idx <- (crec[i]+1):crec[i+1]
                nse[idx,] <- matrix(unlist(data[-(1:4)]), ncol=ncol)
            }, error=function(err) {
                .throw(SRError("Input/Output",
                               sprintf("parsing 'nsePattern'\n  file: %s\n  error: %s",
                                       fls[[i]], conditionMessage(err))))
            })
        }
    }
    readInfo <- ShortReadIntensityInfo(lane=df[[1]], tile=df[[2]],
                                       x=df[[3]], y=df[[4]])
    if (!withNse)
        ShortReadIntensity(intensity=int, readInfo=readInfo)
    else
        ShortReadIntensity(intensity=int, nse=nse, readInfo=readInfo)
}

setMethod(readIntensities, "character", .readIntensities_character)
