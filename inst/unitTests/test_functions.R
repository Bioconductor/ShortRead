## readFastq

test_readFastq_autoDetectType <- function() 
{
    src <- system.file("unitTests","cases", package="ShortRead")
    srq <- readFastq(file.path(src, "sanger.fastq"))
    checkTrue(class(quality(srq)) == "FastqQuality")
    srq <- readFastq(file.path(src, "solexa.fastq"))
    checkTrue(class(quality(srq)) == "SFastqQuality")
    srq <- readFastq(file.path(src, "solexa.fastq"),
                     qualityType="FastqQuality")
    checkTrue(class(quality(srq)) == "FastqQuality")
}

test_readFastq_withids <- function() {
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
    rfq1 <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt",
                      withIds=FALSE)
    checkIdentical(as.character(sread(rfq)), as.character(sread(rfq1)))
    checkIdentical(as.character(quality(quality(rfq))),
                   as.character(quality(quality(rfq1))))
    checkIdentical(as.character(id(rfq1)), character(length(rfq1)))
}

test_readFastq_zerowidth <- function() {
    fl <- tempfile();
    writeLines("@ \n\n+\n", fl)
    fq <- readFastq(fl)
    checkTrue(validObject(fq))
    checkIdentical(0L, width(fq))
}

## alphabetByCycle

checkAlphabetByCycle <- function(obj) {
    abc <- alphabetByCycle(obj)
    validObject(abc)
    checkEquals(length(obj)*unique(width(obj)), sum(abc))
}

test_alphabetByCycle <- function() {
    sp <- SolexaPath(system.file('extdata', package="ShortRead"))
    sq <- readFastq(sp)

    checkAlphabetByCycle(sread(sq))
    checkAlphabetByCycle(quality(sq))

    obj <- alphabetByCycle(sq)
    validObject(obj)
    checkEquals(c(18, 94, 36), dim(obj))

    checkEqualsNumeric(alphabetByCycle(sread(sq)),
                       apply(obj, c(1, 3), sum))
    checkEqualsNumeric(alphabetByCycle(quality(sq)),
                       apply(obj, 2:3, sum))

    obj <- rowSums(alphabetByCycle(id(sq)))
    obj <- obj[obj != 0]
    exp <- table(unlist(strsplit(as.character(id(sq)), ""),
                        use.names=FALSE))
    checkTrue(setequal(names(obj), names(exp)))
    checkIdentical(as.numeric(exp[names(obj)]),
                   as.vector(obj))

    srq <- ShortReadQ(DNAStringSet(), FastqQuality())
    abc <- alphabetByCycle(srq)
    alf <- alphabet(sread(srq))
    qalf <- alphabet(quality(srq))
    checkIdentical(matrix(0L, nrow=length(alf), ncol=0,
                          dimnames=list(alphabet=alf,
                            cycle=character(0))),
                   alphabetByCycle(sread(srq)))
    checkIdentical(array(0L, dim=c(18, 94, 0),
                         dimnames=list(base=alf, quality=qalf,
                           cycle=character(0))),
                   alphabetByCycle(srq))
}

## countLines

test_countLines <- function() {
    sp <- SolexaPath(system.file('extdata', package="ShortRead"))
    nlines <- countLines(analysisPath(sp), "s_1_sequence.txt")
    exp <- 1024; names(exp) <- "s_1_sequence.txt"
    checkEquals(exp, nlines)
    dir <- tempfile()
    dir.create(dir)
    checkException(countLines(dir), silent=TRUE)
}

## sort / order

test_order_stats <- function()
{
    checkIdentical(integer(0), srrank(AlignedRead()))
    checkIdentical(integer(0), srorder(AlignedRead()))
    checkIdentical(logical(0), srduplicated(AlignedRead()))
}

test_alphabetOrder <- function() {
    ## setup
    oldc <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", oldc))
    Sys.setlocale("LC_COLLATE", "C")
    sp <- SolexaPath(system.file('extdata', package='ShortRead'))
    rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")

    checkEquals(srorder(sread(rfq)),
                order(as.character(sread(rfq))))
    checkEquals(srorder(quality(rfq)),
                order(as.character(quality(quality(rfq)))))

    checkEquals(srduplicated(sread(rfq)),
                duplicated(as.character(sread(rfq))))
    checkEquals(srduplicated(quality(rfq)),
                duplicated(as.character(quality(quality(rfq)))))
}

## _mark_field (C code)

test_mark_field <- function() {
    fl <- tempfile()
    do <- function(s, fl) {
        doexp(s, strsplit(unlist(strsplit(s, "\n")), "\t"), fl)
    }
    doexp <- function(s, exp, fl) {
        writeChar(s, fl)
        res <- .Call("_mark_field_test", fl, "\t", c(2L, 3L),
                     PACKAGE="ShortRead")
        checkIdentical(exp, res)
    }

    do("a\tb\tc\nd\te\tf\n", fl)
    do("a\t\tc\nd\te\tf\n", fl)
    do("\tb\tc\nd\te\tf\n", fl)
    do("\t\tc\nd\te\tf\n", fl)

    ## trailing \t are problematic for strsplit
    doexp("a\tb\t\nd\te\tf\n",
          list(c("a","b",""), c("d","e","f")),
          fl)
    doexp("a\t\t\nd\te\tf\n",
          list(c("a","",""), c("d","e","f")),
          fl)

    writeChar("\n", fl)
    res <- .Call("_mark_field_test", fl, "\t", c(1L,1L),
                 PACKAGE="ShortRead")
    checkIdentical(list(""), res)
}
