test_readXStringColumns_toIUPAC <- function() {
    src <- system.file("unitTests", "cases", package="ShortRead") 
    fl <- file.path(src, "s_2_export_toIUPAC.txt")

    colClasses <- rep(list(NULL), 22)
    colClasses[9:10] <- c(sread="DNAString", quality="BString")
    names(colClasses)[9:10] <- c("sread", "quality")
    res <- readXStringColumns(dirname(fl), basename(fl),
                              colClasses=colClasses)

    ## '.' converted to "-" in DNAString, but not BString
    checkTrue(all(gregexpr("-", as.character(res$sread[1]))[[1]] ==
                  c(5, 10)))
    checkTrue(all(gregexpr("\\.", as.character(res$quality[1]))[[1]] == 5))
    checkTrue(all(gregexpr("-", as.character(res$quality[1]))[[1]] == 10))
}

test_readXStringColumns_skip_nrows <- function()
{
    what <- vector("list", 22)
    what[[2]] <- character()
    colClasses <- what
    colClasses[[2]] <- "DNAString"

    ## single file
    dir <- system.file("unitTests", "cases", package="ShortRead"); fl <- "s_1_results_head.txt"
    check <- function(dir, fl, skip, nrows) {
        pth <- file.path(dir, fl)
        exp <- DNAStringSet(scan(pth, what, nmax=nrows, skip=skip,
                                 fill=TRUE, quiet=TRUE)[[2]])
        obs <- readXStringColumns(dir, fl, colClasses=colClasses,
                                  nrows=nrows, skip=skip)[[1]]
        checkEquals(as.character(exp),
                    as.character(obs))
    }
    check(dir, fl, 0L,-1L)
    check(dir, fl, 100L, -1L)
    check(dir, fl, 0L, 100L)
    check(dir, fl, 100L, 100L)

    ## multiple files
    dir <- system.file("unitTests", "cases", package="ShortRead"); pattern <- "s_1_results_head.*txt"
    mcheck <- function(dir, pattern, skip=0L, nrows=-1L) {
        fls <- list.files(dir, pattern, full=TRUE)
        exp <- vector("list",length(fls))
        nread <- 0
        for (i in seq_along(fls)) {
            if (nrows > 0 && nread >= nrows)
                break
            exp[[i]] <- scan(fls[i], what=what, fill=TRUE, skip=skip,
                             nmax=nrows-nread, quiet=TRUE)
            nread <- nread + length(exp[[i]][[2]])
        }
        exp <- DNAStringSet(unlist(exp))
        obs <- readXStringColumns(dir, pattern, colClasses=colClasses,
                                  skip=skip, nrows=nrows)[[1]]
        checkTrue(validObject(obs))
        checkEquals(as.character(exp), as.character(obs))
    }
    mcheck(dir, pattern, 0L)
    mcheck(dir, pattern, 100L)
    mcheck(dir, pattern, 0L, 500L)
    mcheck(dir, pattern, 0L, 1500L)
    mcheck(dir, pattern, 0L, 15000L)
    mcheck(dir, pattern, 100L, 500L)
    mcheck(dir, pattern, 100L, 1500L)
    mcheck(dir, pattern, 100L, 15000L)
}
