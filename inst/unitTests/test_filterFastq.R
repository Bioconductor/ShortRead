sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")

.all_equal <-
    function(target, current, ...)
{
    ac <- as.character
    all.equal(ac(sread(target)), ac(sread(current))) &&
        all.equal(ac(quality(quality(target))),
                  ac(quality(quality(current)))) &&
            all.equal(ac(id(target)), ac(id(current)))
}                         

test_filterFastq <- function() {
    tf <- c(TRUE, FALSE)
    exp <- readFastq(fl)[tf]

    filt <- function(x) x[tf]
    dest <- filterFastq(fl, tempfile(), filter=filt)
    checkTrue(.all_equal(exp, readFastq(dest)))

    dest <- filterFastq(fl, tempfile(), filter=filt, yieldSize=100)
    checkTrue(.all_equal(exp,readFastq(dest)))

    filt <- function(x) tf
    rule <- FilterRules(list(filt=filt))
    dest <- filterFastq(fl, tempfile(), filter=rule)
    checkTrue(.all_equal(exp,readFastq(dest)))

    dest <- tempfile()
    file.create(dest)
    obs <- tryCatch(filterFastq(fl, dest, filter=filt),
                    error=conditionMessage)
    checkIdentical(sprintf("'destinations' exist:\n  %s", dest), obs)
}
