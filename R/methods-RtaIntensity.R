## RtaIntensity

RtaIntensity <-
    function(intensity=array(0, c(0, 0, 0)),
             measurementError=array(0, c(0, 0, 0)),
             readInfo=SolexaIntensityInfo(
               lane=integer()[seq_len(nrow(intensity))]),
             ...)
{
    .hasMeasurementError <- mkScalar(!missing(measurementError))
    new("RtaIntensity",
        intensity=ArrayIntensity(intensity),
        measurementError=ArrayIntensity(measurementError),
        readInfo=readInfo,
        .hasMeasurementError=.hasMeasurementError,
        ...)
}

.readIntensities_RtaIntensity <-
    function(dirPath, pattern=character(0), ...,
             lane=integer(0), cycles=integer(0), cycleIteration=1L,
             tiles=integer(0),
             laneName=sprintf("L%.3d", lane),
             cycleNames=sprintf("C%d.%d", cycles, cycleIteration),
             tileNames=sprintf("s_%d_%d", lane, tiles),
             posNames=sprintf("s_%d_%.4d_pos.txt", lane, tiles),
             withVariability=TRUE, verbose=FALSE)
{
    .check_type_and_length(dirPath, "character", 1)
    .check_type_and_length(pattern, "character", NA)
    .check_type_and_length(lane, "integer", 1)
    .check_type_and_length(cycles, "integer", NA)
    .check_type_and_length(cycleIteration, "integer", 1)
    .check_type_and_length(tiles, "integer", NA)
    
    posFilenames <- file.path(dirPath, posNames)
    ok <- sapply(posFilenames, file.exists)
    if (!all(ok)) {
        msg <-
            sprintf("%d pos files do not exist\n  %s", sum(!ok),
                    paste(selectSome(posFilenames[!ok]), collapse="\n  "))
        .throw(SRError("UserArgumentMismatch", msg))
    }
    if (verbose)
        message("reading 'pos' files")
    readInfo <- do.call(rbind, mapply(function(fl, lane, tile) {
        cbind(lane=lane, tile=tile, read.table(fl, col.names=c("x", "y")))
    }, posFilenames, laneName, tileNames, SIMPLIFY=FALSE, USE.NAMES=FALSE))
    readInfo <- do.call(SolexaIntensityInfo, readInfo)

    laneDirname <- file.path(dirPath, laneName)
    if (!file.exists(laneDirname)) {
        msg <-
            sprintf("unknown lane directory\n  %s", laneDirname)
        .throw(SRError("UserArgumentMismatch", msg))
    }

    cycleDirnames <- file.path(laneDirname, cycleNames)
    ok <- sapply(cycleDirnames, file.exists)
    if (!all(ok)) {
        msg <-
            sprintf("%d cycle directories do not exist\n  %s", sum(!ok),
                    paste(selectSome(cycleDirnames[!ok]), collapse="\n  "))
        .throw(SRError("UserArgumentMismatch", msg))
    }

    if (verbose)
        message("reading 'cif' files")
    cif <- .read_cif_or_cnf(cycleDirnames, tileNames, ".cif")
    if (withVariability) {
        if (verbose)
            message("reading 'cnf' files")
        cnf <- .read_cif_or_cnf(cycleDirnames, tileNames, ".cnf")
        RtaIntensity(intensity=cif, measurementError=cnf,
                     readInfo=readInfo)
    } else {
        RtaIntensity(intensity=cif, readInfo=readInfo)
    }
}

.read_cif_or_cnf_file <-
    function(fileName)
{
    conn <- file(fileName, "rb")
    on.exit(close(conn))

    ## header
    id <- rawToChar(readBin(conn, "raw", 3L))
    if (id != "CIF")
        stop("not a CIF / CNF file:\n  id: ", id,
             "\n  file: ", fileName)
    version <- readBin(conn, "integer", 1L, 1L, signed=FALSE)
    if (version != 1L)
        stop("unknown CIF / CNF version:\n  version: ", version,
             "\n  file: ", fileName)
    dataType <- readBin(conn, "integer", 1L, 1L, signed=FALSE)
    firstCycle <- readBin(conn, "integer", 1L, 2L,
                          signed=FALSE, endian="little")
    numberOfCycles <- readBin(conn, "integer", 1L, 2L,
                              signed=FALSE, endian="little")
    numberOfClusters <- readBin(conn, "integer", 1L, 4L,
                                signed=FALSE, endian="little")

    ## data
    m <- readBin(conn, "integer", 4 * numberOfClusters,
                 dataType, endian="little")
    if (length(m) != 4 * numberOfClusters)
        stop("incorrect number of CIF data values:",
             "\n  expected: ", 4 * numberOfClusters,
             "\n  found: ", length(m),
             "\n  file: ", fileName)
    m
}

.read_cif_or_cnf <-
    function(cycleDirs, tileNamesRoot, ext)
{
    res <- lapply(cycleDirs, function(dir, tileNames) {
        fls <- file.path(dir, tileNames)
        unlist(lapply(fls, .read_cif_or_cnf_file))
    }, paste(tileNamesRoot, ext, sep=""))
    isNull <- sapply(res, is.null)
    if (any(isNull))
        stop("no CIF or CNF files matching pattern",
             "\n  pattern: '", paste(tileNamesRoot, ext, sep=""), "'",
             "\n  directories:\n    ",
             paste(cycleDirs[isNull], collapse="\n    ", sep=""))
    nClusters <- unique(sapply(res, length)) / 4
    if (length(nClusters) != 1L)
        stop("cluster counts differ between cycles",
             "\n  found: ", paste(nClusters, collapse=" "))
    nms <- list(NULL, c("A", "C", "G", "T"), basename(cycleDirs))
    array(unlist(res), dim=c(nClusters, 4L, length(cycleDirs)),
          dimnames=nms)
}
