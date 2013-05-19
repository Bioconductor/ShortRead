sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")
rfq <- readFastq(fl)

.check <- function(xexp, xobs)
    checkIdentical(as.character(xexp), as.character(xobs))

test_trimTails_BStringSet <- function()
{
    .check(BStringSet("CCCBBB"),
                     trimTails(BStringSet("CCCBBBAAA"), 1, "A"))
    .check(BStringSet("CCCABBB"),
                     trimTails(BStringSet("CCCABBBAAA"), 2, "A"))
    .check(BStringSet("CCCABBBAB"),
           trimTails(BStringSet("CCCABBBABAA"), 2, "A", successive=TRUE))
    .check(BStringSet("CCC"),
           trimTails(BStringSet("CCCABBBABAA"), 2, "B", successive=TRUE))
    .check(BStringSet(), trimTails(BStringSet("CCCABBBABAA"), 1, "C"))
}

test_trimTails_QualityScore <- function()
{
    .qq <- function(x) quality(quality(x))
    checkTrue(validObject(trimTails(rfq, 1, "H")))
    .check(.qq(rfq), .qq(trimTails(rfq, 1, " ")))
    .check(BStringSet(), .qq(trimTails(rfq, 1, "]")))
}

test_trimTails_XStringQuality <- function()
{
    .qq <- function(x) quality(quality(x))
    .qb <- function(x) as(x, "BStringSet")
    qual <- as(quality(rfq), "PhredQuality")
    checkTrue(validObject(trimTails(qual, 1, "H")))
    .check(.qq(rfq), .qb(trimTails(qual, 1, "!")))
    .check(BStringSet(), .qb(trimTails(qual, 1, "]")))
}

test_trimTails_file <- function()
{
    exp <- width(trimTails(rfq, 1, "H"))
    dest <- trimTails(fl, 1, "H", destinations=tempfile())
    checkIdentical(exp, width(readFastq(dest)))
}

test_trimTailw <- function()
{
    b <- BStringSet("BBBBBB")
    checkIdentical(BStringSet(), trimTailw(b, 1L, "C", 3L))
    checkIdentical(BStringSet(), trimTailw(b, 1L, "B", 3L))
    checkIdentical(b, trimTailw(b, 1L, "A", 3L))
    checkIdentical(BStringSet(), trimTailw(b, 3L, "C", 1L))
    checkIdentical(b, trimTailw(b, 4L, "C", 1L))

    b <- BStringSet("DDDBBBBB")
    checkIdentical(BStringSet("DDD"), trimTailw(b, 2L, "C", 1L))
    checkIdentical(BStringSet("DD"), trimTailw(b, 1L, "C", 1L))
    checkIdentical(BStringSet("D"), trimTailw(b, 1L, "C", 2L))
    checkIdentical(BStringSet(), trimTailw(b, 1L, "C", 3L))

    b <- BStringSet("DDDBDBBBB")
    checkIdentical(BStringSet("DDDBD"), trimTailw(b, 2L, "C", 1L))
    checkIdentical(BStringSet("DDDBDB"), trimTailw(b, 4L, "C", 2L))
}

test_trimTailw_file <- function()
{
    exp <- width(trimTailw(rfq, 2L, "C", 1L))
    dest <- trimTailw(fl, 2L, "C", 1L, destinations=tempfile())
    checkIdentical(exp, width(readFastq(dest)))
}
