sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")
rfq <- readFastq(fl)

test_trimEnds <- function()
{
    exp <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 10, 16, 72, 152)
    checkIdentical(as.integer(exp), tabulate(width(trimEnds(rfq, "I"))))

    rng <- trimEnds(sread(rfq), "G", relation="==", ranges=TRUE)
    checkTrue(!all(1L == start(rng)))
    checkTrue(!all(end(rfq) == end(rng)))

    checkTrue(all(1L == start(trimEnds(sread(rfq), "G", left=FALSE,
                    relation="==", ranges=TRUE))))
    checkTrue(all(width(rfq) == end(trimEnds(sread(rfq), "G", right=FALSE,
                         relation="==", ranges=TRUE))))

    exp <- c(1L, 1L, 3L, 3L, 8L, 8L, 12L, 10L, 41L, 41L, 38L, 43L,
             47L)
    obs <- trimEnds(sread(rfq), c("G", "T"), relation="==")
    checkIdentical(exp, as.vector(table(width(obs))))
}

test_trimEnds_unknown_a <- function()
{
    checkIdentical(as.character(sread(rfq)),
                   suppressWarnings(as.character(trimEnds(sread(rfq), "Z"))))
    obs <- tryCatch(trimEnds(sread(rfq), "Z"), warning=conditionMessage)
    checkIdentical("some 'a' not in alphabet(object)", obs)
}

test_trimEnds_classes <- function()
{
    rng <- trimEnds(quality(rfq), "I", ranges=TRUE)
    checkIdentical(as.character(quality(narrow(quality(rfq),
                                               start(rng), end(rng)))),
                   as.character(quality(trimEnds(quality(rfq), "I"))))
    ## FIXME: additional, e.g., PhredQuality
}

test_trimEnds_file <- function() {
    dest <- trimEnds(fl, "I", destinations=tempfile())
    checkIdentical(width(trimEnds(rfq, "I")), width(readFastq(dest)))
}
