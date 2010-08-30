test_SRList_construction <- function() {
    srl <- SRList()
    checkTrue(validObject(srl))
    checkEquals(0, length(srl))

    srl <- SRList(list())
    checkTrue(validObject(srl))
    checkEquals(0, length(srl))

    srl <- SRList(list(1))
    checkTrue(validObject(srl))
    checkEquals(1, length(srl))

    srl <- SRList(list(1, 2))
    checkTrue(validObject(srl))
    checkEquals(2, length(srl))
    checkEquals(1, length(srl[[1]]))

    srl <- SRList(1, 2)
    checkTrue(validObject(srl))
    checkEquals(2, length(srl))
    checkEquals(1, length(srl[[1]]))
}
