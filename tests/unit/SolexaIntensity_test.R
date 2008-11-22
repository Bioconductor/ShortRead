test_SolexaIntensity_construction <- function()
{
    checkTrue(validObject(SolexaIntensityInfo()))
    checkTrue(validObject(SolexaIntensityInfo(lane=rep(1, 10))))

    checkTrue(validObject(SolexaIntensity()))
    checkTrue(validObject(SolexaIntensity(intensity=array(0,c(1,2,3)))))
    checkException(SolexaIntensity(measurementError=array(0,c(1,2,3))),
                   silent=TRUE)

    checkTrue(validObject(SolexaIntensity(intensity=array(0,c(1,2,3)),
                                          measurementError=array(0,c(1,2,3)))))

    checkException(SolexaIntensityInfo()[,"lane"],
                   silent=TRUE)
}

test_SolexaIntensity_access <- function()
{
    checkException(measurementError(SolexaIntensity()),
                   silent=TRUE)
}

test_SolexaIntensity_io <- function() 
{
    sp <- SolexaPath(system.file("extdata", package="ShortRead"))
    int <- readIntensities(sp)
    checkIdentical(c(256L, 4L, 36L), dim(intensity(int)))
    checkIdentical(c(256L, 4L, 36L), dim(measurementError(int)))
    checkIdentical(256L, nrow(pData(readInfo(int))))

    int <- readIntensities(sp, withVariability=FALSE)
    checkIdentical(c(256L, 4L, 36L), dim(intensity(int)))
    checkIdentical(256L, nrow(pData(readInfo(int))))
    checkException(measurementError(int), silent=TRUE)
}
