#include <string.h>
#include <R.h>
#include <Rinternals.h>

SEXP pileup( SEXP start, SEXP fraglength, SEXP chrlength, SEXP dir, SEXP readlength,
   SEXP offset )
{
   SEXP res;
   int i, j, st, end, offs;
  
   offs = INTEGER(offset)[0];
   PROTECT( res = allocVector( INTSXP, INTEGER(chrlength)[0] ) );
   memset( INTEGER(res), 0, length(res) * sizeof(int) );

   for( i = 0; i < length(start); i++ )
      if( INTEGER(dir)[ length(dir) == 1 ? 0 : i] == 2 ) {
         /* forward direction */
	 end = INTEGER(start)[i] + INTEGER(fraglength)[ length(fraglength) == 1 ? 0 : i];
	 if( end - offs > length(res) )
	    error( "'chrlength' is too small" );	 
         for( j = INTEGER(start)[i]; j < end; j++ )
	    INTEGER(res)[j-offs] += 1;
      } else {
         /* backward strand */
	 st = INTEGER(start)[i] + INTEGER(readlength)[ length(readlength) == 1 ? 0 : i] - 1;
	 if( st - offs >= length(res) )
	    error( "'chrlength' is too small" );	 
	 end = st - INTEGER(fraglength)[ length(fraglength) == 1 ? 0 : i];
	 if( end - offs < 0 )
	    error( "Lower bound of pile-up vector exceeded." );	 
         for( j = st; j > end; j-- )
	    INTEGER(res)[j-offs] += 1;
      }   
      
   UNPROTECT( 1 );  
   return res;
}
