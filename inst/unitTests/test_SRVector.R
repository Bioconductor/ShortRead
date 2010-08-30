test_SRVector_construction <- function() {
    check <- function(srv, cls, len) {
        checkTrue(validObject(srv))
        checkEquals(cls, vclass(srv))
        checkEquals(len, length(srv))
    }
        
    check(SRVector(vclass="numeric"), "numeric", 0)
    check(SRVector(1), "numeric", 1)

    checkException(SRVector(),silent=TRUE)
    checkException(SRVector("a", vclass="numeric"), silent=TRUE)
    checkException(SRVector(1, "a"), silent=TRUE)
}
