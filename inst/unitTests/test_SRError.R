throw <- ShortRead:::.throw

## SRError

test_SRError_construction <- function() {
    checkTrue(validObject(SRError("UnspecifiedError", "message")))

    checkException(SRError(), silent=TRUE)
    checkException(SBError("UnspecifiedError"),
                   silent=TRUE) # must have message
    checkException(SRError("Bad error class", "Message"),
                   silent=TRUE) # must have valid class
}

test_SRError_throw <- function() {
    err <- SRError("UnspecifiedError", "error message")
    checkException(throw(err), silent=TRUE)
}

test_SRError_throw <- function() {
    err <- SRError("UnspecifiedError", "error message")
    checkTrue(tryCatch(throw(err),
                       SRError=function(err) TRUE))
    checkTrue(tryCatch(throw(err),
                       UnspecifiedError=function(err) TRUE))
}

## SRWarn

test_SRWarn_construction <- function() {
    checkTrue(validObject(SRWarn("UnspecifiedWarning", "message")))

    checkException(SRWarn(), silent=TRUE)
    checkException(SRWarn("UnspecifiedWarning"), silent=TRUE)
    checkException(SRWarn("Bad Warn Class"), silent=TRUE)
}

test_SRWarn_throw <- function() {
    old.opt <- options(warn=2)
    on.exit(options(old.opt))
    warn <- SRWarn("UnspecifiedWarning", "warning message")
    checkException(throw(warn), silent=TRUE)
}

test_SRWarn_catch <- function() {
    old.opt <- options(warn=2)
    on.exit(options(old.opt))
    warn <- SRWarn("UnspecifiedWarning", "warning message")

    checkTrue(tryCatch(throw(warn),
                       SRWarn=function(warn) TRUE))
    checkTrue(tryCatch(throw(warn),
                       UnspecifiedWarning=function(warn) TRUE))
}
