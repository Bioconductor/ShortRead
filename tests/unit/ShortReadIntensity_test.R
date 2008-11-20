sp <- SolexaPath(system.file("extdata", package="ShortRead"))

test_ShortReadIntensity_construction <- function()
{
    checkTrue(validObject(ShortReadIntensityInfo()))
    checkTrue(validObject(ShortReadIntensityInfo(lane=rep(1, 10))))

    checkTrue(validObject(ShortReadIntensity()))
    checkTrue(validObject(ShortReadIntensity(intensity=matrix(0,1,2))))
    checkException(ShortReadIntensity(nse=matrix(0,1,2)),
                   silent=TRUE)

    checkTrue(validObject(ShortReadIntensity(intensity=matrix(0,1,2),
                                             nse=matrix(0,1,2))))

    slim <- ShortReadIntensityInfo()[,"lane"]
    checkException(ShortReadIntensity(readInfo=slim),
                   silent=TRUE)
}

test_ShortReadIntensity_access <- function()
{
    checkException(nse(ShortReadIntensity()),
                   silent=TRUE)
}

test_ShortReadIntensity_io <- function() 
{
    int <- readIntensities(sp)
    checkIdentical(c(256L, 144L), dim(intensity(int)))
    checkIdentical(c(256L, 144L), dim(nse(int)))
    checkIdentical(256L, nrow(pData(readInfo(int))))

    int <- readIntensities(sp, withNse=FALSE)
    checkIdentical(c(256L, 144L), dim(intensity(int)))
    checkIdentical(256L, nrow(pData(readInfo(int))))
    checkException(nse(int), silent=TRUE)
}
