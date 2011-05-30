sp <- SolexaPath(system.file('extdata', package='ShortRead'))
rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")

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
    qual <- as(quality(rfq), "PhredQuality")
    checkTrue(validObject(trimTails(rfq, 1, "H")))
    .check(.qq(rfq), .qq(trimTails(rfq, 1, " ")))
    .check(BStringSet(), .qq(trimTails(rfq, 1, "]")))
}
