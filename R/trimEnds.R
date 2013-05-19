setMethod("trimTails", "character",
    function(object, k, a, successive=FALSE, ..., destinations,
             ranges=FALSE)
{
    filterFastq(object, destinations, k=k, a=a, ..., filter=trimTails,
                ranges=ranges)
})

setMethod("trimTailw", "character",
    function(object, k, a, halfwidth, ..., destinations, ranges=FALSE)
{
    filterFastq(object, destinations, k=k, a=a, halfwidth=halfwidth,
                ..., filter=trimTailw, ranges=ranges)
})

setMethod("trimEnds", "character",
    function(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
                    ..., destinations, ranges=FALSE)
{
    filterFastq(object, destinations, a=a, left=left, right=right,
                relation=relation, ..., filter=trimEnds,
                ranges=ranges)
})
