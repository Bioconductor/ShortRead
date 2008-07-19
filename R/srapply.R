.fapply <- function(...) {
    catchErrs <- function(FUN) {
        function(...) {
            tryCatch({
                FUN(...)
            }, error=function(err) {
                SRError("RemoteError",
                        paste(capture.output(conditionCall(err)),
                              conditionMessage(err), sep="\n  "))
            })
        }
    }
    if (is.loaded("mpi_comm_size")) {
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
                res <- parLapply(X, CFUN, ..., verbose=verbose)
                res
            }
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

.reduce <- function() {
    function(lst, ...) {
        errs <- sapply(lst, is, "SRError")
        if (any(errs)) {
            elts <- lst[errs]
            msg <- paste(sapply(elts, .type), sapply(elts, .message),
                         sep=": ", collapse="\n  ")
            .throw(SRWarn("RemoteWarning",
                          "elements: %s\n  %s",
                          paste(which(errs), collapse=" "),
                          msg))
            lst <- lst[!errs]
        }
        lst
    }
}

srapply <- function(X, FUN, ...,
                    fapply=.fapply(), reduce=.reduce(), verbose=FALSE) {
    result <- fapply(X, FUN, ..., verbose=verbose)
    reduce(result)
}
