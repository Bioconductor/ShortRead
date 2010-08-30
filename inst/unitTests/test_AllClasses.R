test_AllClasses_STRAND_LEVELS <- function()
{
    checkIdentical(ShortRead:::.STRAND_LEVELS, levels(strand()))
}

test_AllClasses_new <- function() {
    nmspace <- getNamespace("ShortRead")
    nms <- names(slot(getClass(".ShortReadBase", where=nmspace),
                      "subclasses"))
    ## 'new' with no additional arguments
    ok <- Map(function(x) validObject(new(x)),
              Filter(function(x) !slot(getClass(x, where=nmspace), "virtual"),
                     nms))
    checkTrue(all(unlist(ok)))
}
