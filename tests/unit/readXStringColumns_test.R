test_readXStringColumns_toIUPAC <- function() {
    fl <- file.path("cases/s_2_export_toIUPAC.txt")

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
