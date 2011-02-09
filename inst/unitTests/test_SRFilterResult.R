mkScalar <- Biobase::mkScalar

test_SRFilterResult_constructor <- function()
{
    checkTrue(validObject(SRFilterResult()))

    fr <- SRFilterResult(TRUE)
    checkTrue(validObject(fr))
    checkIdentical(mkScalar(NA_character_), name(fr))
    df <- data.frame(Name=NA_character_, Input=1L, Passing=1L,
                     Op=NA_character_, stringsAsFactors=FALSE)
    checkIdentical(df, stats(fr))

    fr <- SRFilterResult(c(TRUE, FALSE))
    df <- data.frame(Name=NA_character_, Input=2L, Passing=1L,
                     Op=NA_character_, stringsAsFactors=FALSE)
    checkIdentical(df, stats(fr))

    fr <- SRFilterResult(c(TRUE, FALSE),"A")
    df <- data.frame(Name="A", Input=2L, Passing=1L, Op=NA_character_,
                     stringsAsFactors=FALSE)
    checkIdentical(df, stats(fr))
}

test_SRFilterResult_logic <- function()
{
    a <- SRFilterResult(c(TRUE, FALSE), "A")
    b <- SRFilterResult(c(FALSE, TRUE), "B")

    checkTrue(all(a|b))
    exp <-
        structure(list(Name = c("A", "B", "(A | B)"),
                       Input = c(2L, 2L, 2L), Passing = c(1L, 1L, 2L),
                       Op = c(NA, NA, "|")),
                  .Names = c("Name", "Input", "Passing", "Op"),
                  row.names = c(NA, -3L), class = "data.frame")
    checkIdentical(exp, stats(a|b))

    checkTrue(all(!b == !(b@.Data)))
    exp <-
        structure(list(Name = c("B", "!(B)"), Input = c(2L, 2L),
                       Passing = c(1L, 1L), Op = c(NA, "!")),
                  .Names = c("Name", "Input", "Passing", "Op"),
                  row.names = c(NA, -2L), class = "data.frame")
    checkIdentical(exp, stats(!b))
}

test_SRFilterResult_SRFilter <- function()
{
    fa <- srFilter(function(x) logical(length(x)), "A")
    x <- fa(1:10)
    checkIdentical(logical(10L), x@.Data)
    checkIdentical(mkScalar("A"), name(x))
    exp <- 
        structure(list(Name = "A", Input = 10L, Passing = 0L, Op =
                       NA_character_),
                  .Names = c("Name", "Input", "Passing", "Op"),
                  row.names = c(NA, -1L), class = "data.frame")
    checkIdentical(exp, stats(x))

    x <- fa(1:10) & fa(1:10)
    checkIdentical(logical(10L), x@.Data)
    checkIdentical(mkScalar("(A & A)"), name(x))
    exp <-
        structure(list(Name = c("A", "A", "(A & A)"),
                       Input = c(10L, 10L, 10L), Passing = c(0L, 0L, 0L),
                       Op = c(NA, NA, "&")),
                  .Names = c("Name", "Input", "Passing", "Op"),
                  row.names = c(NA, -3L), class = "data.frame")
    checkIdentical(exp, stats(x))

    fb <- srFilter(function(x) !logical(length(x)), "B")
    x <- fa(1:10) | fb(1:10)
    checkIdentical(!logical(10L), x@.Data)
    checkIdentical(mkScalar("(A | B)"), name(x))
    exp <-
        structure(list(Name = c("A", "B", "(A | B)"),
                       Input = c(10L, 10L, 10L), Passing = c(0L, 10L, 10L),
                       Op = c(NA, NA, "|")),
                  .Names = c("Name", "Input", "Passing", "Op"),
                  row.names = c(NA, -3L), class = "data.frame")
    checkIdentical(exp, stats(x))
}
