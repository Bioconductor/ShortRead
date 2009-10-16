pileup <- function (
   start, fraglength, chrlength, dir = strand( "+" ),
   readlength = fraglength,
   offset = 1 )
{
   .Deprecated("coverage")
   stopifnot( is.factor(dir) )
   stopifnot( length(levels(dir)) == 3 &&
             all( levels(dir) == c("-", "+", "*") ) )
   stopifnot( is.numeric(start) &&
             is.numeric(fraglength) && is.numeric(chrlength) 
      && is.numeric(readlength) )
   stopifnot( all(!is.na(start)) && all(!is.na(fraglength)) &&
              all(!is.na(chrlength)) && all(!is.na(readlength)) )
   stopifnot( length(fraglength) == length(start) || length(fraglength) == 1 )
   stopifnot( length(chrlength) == 1 )
   stopifnot( length(readlength) == length(start) || length(readlength) == 1 )
   stopifnot( length(dir) == length(start) || length(dir) == 1 )
   stopifnot( length(offset) == 1 )
   	    
   .Call(.pileup, as.integer(start), as.integer(fraglength), as.integer(chrlength),
       as.integer(dir), as.integer(readlength), as.integer(offset) )
}
