sp <- SolexaPath(system.file("extdata", package="ShortRead"))
sr <- as(readFastq(sp, "s_1_sequence.txt"), "ShortRead")
    
.equals <- function(x, y)
{
    checkIdentical(as.character(sread(x)), as.character(sread(y)))
    checkIdentical(as.character(id(x)), as.character(id(y)))
}

test_ShortRead_construction <- function() {
    obj <- ShortRead()
    checkTrue(class(obj) == "ShortRead")
    checkTrue(validObject(obj))

    obj <- ShortRead(sread(sr))
    checkTrue(class(obj) == "ShortRead")
    checkTrue(validObject(obj))
    .equals(new("ShortRead", sread=DNAStringSet(sread(sr)),
                id=BStringSet(rep("", length(sr)))), obj)

    obj <- ShortRead(sread(sr), id(sr))
    checkTrue(class(obj) == "ShortRead")
    checkTrue(validObject(obj))
    .equals(sr, obj)
}

test_ShortRead_narrow <- function() {
    obj <- narrow(sr, start=1, end=10)
    checkTrue(class(obj) == "ShortRead")
    checkTrue(length(obj) == length(sr))
    checkTrue(unique(width(obj)) == 10)
    checkIdentical(as.character(sread(obj)),
                   substr(as.character(sread(sr)), 1, 10))

    checkIdentical(narrow(sr, start=start(sread(sr))), sr)
}
