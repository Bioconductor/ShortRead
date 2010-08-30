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

test_IparIntensity_io <- function()
{
    src <- system.file("unitTests","cases",package="ShortRead") 
    int <- readIntensities(src, type="IparIntensity",
                           intExtension="_int_head.txt.p",
                           nseExtension="_nse_head.txt.p",
                           posExtension="_pos_head.txt")
    checkIdentical(c(5L, 4L, 3L), dim(int))

    checkIdentical(structure(c(11.7, 49.4, 14.3, 110.6, -1.1, 3.2,
                               86.2, 5.9, 18.8, 120.2, 21.1, 218.8,
                               96.2, 2.7, 177.4, 9.1, 55.9, 340.9,
                               112, 164.8, 0.8, 2.6, 15.7, 4.6), .Dim
                               = c(2L, 4L, 3L), .Dimnames = list(
                               NULL, c("A", "C", "G", "T"), NULL)),
                   as(intensity(int), "array")[1:2,,])

    checkIdentical(structure(c(8.5, 6.2, 9, 9.2, 3, 3.1, 6.6, 6.4,
                               12.6, 11.9, 12.5, 11, 3.8, 3.6, 6.6,
                               5.8, 13.6, 10.5, 12.4, 12.7, 4.7, 4.2,
                               7.3, 6.4), .Dim = c(2L, 4L, 3L),
                               .Dimnames = list(NULL, c("A", "C", "G",
                               "T"), NULL)),
                   as(measurementError(int), "array")[1:2,,])

    checkIdentical(structure(list(lane = structure(c(1L, 1L, 1L, 1L, 1L),
                                    class = "factor", .Label = "1"), 
                                  tile = c(1L, 1L, 1L, 1L, 1L),
                                  x = c(-0.47, -0.45, -0.45, -0.44,
                                    -0.43),
                                  y = c(1073.78, 1558.67, 1157.37,
                                    144.35, 1497.99 )),
                             .Names = c("lane", "tile", "x", "y"),
                             row.names = c(NA, -5L), class = "data.frame"),
                   pData(readInfo(int)))
}
