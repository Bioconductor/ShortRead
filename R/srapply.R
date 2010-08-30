.fapply <- function(...) {
    catchErrs <- function(FUN) {
        function(...) {
            tryCatch({
                FUN(...)
            }, error=function(err) {
                msg <- paste(capture.output(conditionCall(err)),
                             conditionMessage(err), sep="\n ")
                if (is.loaded("mpi_comm_size"))
                    SRError("RemoteError", msg)
                else
                    SRError("UnspecifiedError", msg)
            })
        }
    }
    if (is.loaded("mpi_comm_size", PACKAGE="Rmpi")) {
        ## 'get()' are to quieten R CMD check, and for no other reason
        commSize <- get("mpi.comm.size", mode="function")
        remoteExec <- get("mpi.remote.exec", mode="function")
        bcastRobj <- get("mpi.bcast.Robj", mode="function")
        parLapply <- get("mpi.parLapply", mode="function")
        function(X, FUN, ..., verbose=FALSE) {
            CFUN <- catchErrs(FUN)
            if (commSize()==0) {
                if (verbose)
                    message("Rmpi loaded, but mpi.comm.size==0; using lapply")
                lapply(X, CFUN, ..., verbose=verbose)
            } else {
                if (verbose)
                    message("using mpi.parLapply")
                libOk <- remoteExec(require, "ShortRead")
                if (!all(unlist(libOk)))
                    .throw(SRError("RemoteError",
                                   "could not 'require(ShortRead)' on %s ",
                                   paste(names(libOk)[!unlist(libOk)],
                                         collapse=", ")))
                wd <- getwd()
                remoteExec(setwd, wd, ret=FALSE)
                if (identical(globalenv(), environment(FUN)))
                    bcastRobj(FUN)
                parLapply(X, CFUN, ..., verbose=verbose)
            }
        }
    } else if (is.loaded("mc_fork", PACKAGE="multicore")) {
        mcLapply <- get('mclapply', envir=getNamespace('multicore'))
        function(X, FUN, ..., verbose=FALSE) {
            CFUN <- catchErrs(FUN)
            if (verbose)
                message("using 'mclapply'")
            mcLapply(X, CFUN, ..., verbose=verbose)
        }
    } else {
        function(X, FUN, ..., verbose=FALSE) {
            CFUN <- catchErrs(FUN)
            if (verbose)
                message("using lapply")
            lapply(X, CFUN, ..., verbose=verbose)
        }
    }
}

.reduce <- function(.minimum_length=0L) {
    function(lst, ...) {
        errs <- sapply(lst, is, "SRError")
        if (any(errs)) {
            elts <- lst[errs]
            msg <- paste(sapply(elts, .type), sapply(elts, .message),
                         sep=": ", collapse="\n  ")
            type <- 
                if (is.loaded("mpi_comm_size", PACKAGE="Rmpi") |
                    is.loaded("mc_fork", PACKAGE="multicore"))
                {
                    "RemoteWarning"
                } else "UnspecifiedWarning"
            .throw(SRWarn(type,
                          "elements: %s\n  %s",
                          paste(which(errs), collapse=" "),
                          msg))
            lst <- lst[!errs]
        }
        if (.minimum_length > length(lst))
            .throw(SRError("ValueUnavailable",
                           "%d elements returned; expected >=%d",
                           length(lst), .minimum_length))
        lst
    }
}

srapply <- function(X, FUN, ...,
                    fapply=.fapply(), reduce=.reduce(),
                    USE.NAMES=FALSE, verbose=FALSE) {
    result <- fapply(X, FUN, ..., verbose=verbose)
    if (USE.NAMES && is.character(X) && is.null(names(result)))
        names(result) <- X
    reduce(result)
}
