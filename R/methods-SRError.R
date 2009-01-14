.SRErrorWarning_types <- c("SRVectorClassDisagreement",
                           "Input/Output",
                           "UserSubset",
                           "UserArgumentMismatch")

.SRError_types <- c("UnspecifiedError",
                    "InternalError",
                    "RemoteError",
                    "InvalidReadFilter",
                    "IncompatibleTypes",
                    "ValueUnavailable",
                    .SRErrorWarning_types)

.SRWarn_types <- c("UnspecifiedWarning",
                   "RemoteWarning",
                   .SRErrorWarning_types)

## Error

setMethod(".srValidity", "SRError", function(object) {
    msg <- NULL
    type <- .type(object)
    if (!type %in% .SRError_types)
        msg <- c(msg, sprintf("'%s' must be one of '%s'",
                              '.type',
                              paste(.SRError_types, collapse="' '")))
    if (is.null(msg)) TRUE else msg
})

SRError <- function(type, fmt, ...) {
    new("SRError", .type=type, .message=sprintf(fmt, ...))
}

.make_getter(slotNames("SRError"))

setMethod(".throw", "SRError", function(object, call=NULL, ...) {
    class <- c(.type(object), "SRError", "error", "condition")
    msg <- paste(.type(object), .message(object), sep="\n  ")
    cond <- structure(list(message=msg, call=call), class=class)
    stop(cond)
})

## Warning

setMethod(".srValidity", "SRWarn", function(object) {
    msg <- NULL
    type <- .type(object)
    if (!type %in% .SRWarn_types)
        msg <- c(msg, sprintf("'%s' must be one of '%s'",
                              '.type', .SRWarn_types))
    if (is.null(msg)) TRUE else msg
})

SRWarn <- function(type, fmt, ...) {
    new("SRWarn", .type=type, .message=sprintf(fmt, ...))
}

setMethod(".throw", "SRWarn", function(object, call=NULL, ...) {
    class <- c(.type(object), "SRWarn", "warning", "condition")
    msg <- paste(.type(object), .message(object), sep="\n  ")
    cond <- structure(list(message=msg, call=call), class=class)
    warning(cond)
})
