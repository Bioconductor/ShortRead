test_IntensityMeasure_subset <- function()
{
    a <- array(1:1000, c(10, 10, 10))
    x <- ArrayIntensity(a)
    checkTrue(all(a==x))

    checkTrue(all(a[,,]==x[,,]))
    checkTrue(all(a[1:5,,]==x[1:5,,]))
    checkTrue(all(a[,1:5,]==x[,1:5,]))
    checkTrue(all(a[1:5,1:5,]==x[1:5,1:5,]))

    checkTrue(all(a[,,1:5]==x[,,1:5]))
    checkTrue(all(a[1:5,,1:5]==x[1:5,,1:5]))
    checkTrue(all(a[,1:5,1:5]==x[,1:5,1:5]))
    checkTrue(all(a[1:5,1:5,1:5]==x[1:5,1:5,1:5]))
}

test_Intensity_subset <- function()
{
    check <- function(obj, m, adf)
    {
        checkTrue(all(m==intensity(obj)))
        checkTrue(all(m==measurementError(obj)))
        checkIdentical(adf, readInfo(obj))
    }

    m <- array(1:1000, c(10, 10, 10))
    adf <- SolexaIntensityInfo(1:10)
    si <- SolexaIntensity(intensity=m, measurementError=m, readInfo=adf)

    ridx <- sample(nrow(m), 5)
    cidx <- sample(ncol(m), 5)
    x <- ArrayIntensity(m)
    checkTrue(all(intensity(si)==x))
    checkTrue(all(measurementError(si)==x))
    check(si[,,], m[,,], adf)
    check(si[ridx,,], m[ridx,,], adf[ridx,])
    check(si[,cidx,], m[,cidx,], adf)
    check(si[ridx, cidx,], m[ridx, cidx,], adf[ridx,])
}
