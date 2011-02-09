SRFilterResult <-
    function(x=logical(), name=NA_character_, input=length(x),
             passing=sum(x), op=NA_character_)
{
    new("SRFilterResult", x, name=mkScalar(as.character(name)[length(name)]),
        stats=data.frame(Name=as.character(name), Input=input, Passing=passing,
          Op=op, stringsAsFactors=FALSE))
}

setMethod(name, "SRFilterResult", function(x, ...) slot(x, "name"))

setMethod(stats, "SRFilterResult", function(x, ...) slot(x, "stats"))

setMethod("Logic", c("SRFilterResult", "SRFilterResult"),
    function(e1, e2)
{
    x <- callNextMethod()
    s1 <- stats(e1); s2 <- stats(e2)
    op <- as.character(.Generic)
    name <- sprintf("(%s %s %s)", name(e1), op, name(e2))
    s <- rbind(stats(e1), stats(e2),
               data.frame(Name=name, Input=length(x), Passing=sum(x),
                          Op=op, stringsAsFactors=FALSE))
    SRFilterResult(x, s$Name, s$Input, s$Passing, s$Op)
})

setMethod("!", "SRFilterResult",
    function(x)
{
    name <- sprintf("!(%s)", name(x))
    y <- callNextMethod()
    s <- rbind(stats(x),
               data.frame(Name=name, Input=length(y), Passing=sum(y),
                          Op="!", stringsAsFactors=FALSE))
    SRFilterResult(y, s$Name, s$Input, s$Passing, s$Op)
})

setMethod(show, "SRFilterResult",
    function(object) 
{
    cat("class:", class(object), "\n")
    cat("name:", name(object), "\n")
    cat("output:", selectSome(object), "\n")
    cat("stats:\n")
    print(stats(object))
})

