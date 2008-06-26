#include <stdio.h>
#include <deque>
#include <string>
#include <R.h>
#include <Rinternals.h>
#include "maqmap.h"

struct seq_meta_info {
   seq_meta_info( int len_, char * name_ )
     : len(len_), name(name_) {};
   int len;
   std::string name;
};

extern "C" SEXP readBfaToc( SEXP bfa_filename )
{
   FILE * fp;
   int name_len, seq_ori_len, seq_len;
   char seq_name[201];
   std::deque< seq_meta_info > seqs;
   if( (! isString(bfa_filename) ) || ( length(bfa_filename) != 1 ) )
      error( "First argument invalid: should be the filename." );

   fp = fopen( CHAR(STRING_ELT(bfa_filename,0)), "r" );
   if( !fp ) {
      char buf[300];
      snprintf( buf, 300, "Failed to open file '%s': %s (errno=%d)",
         CHAR(STRING_ELT(bfa_filename,0)), strerror(errno), errno );
      error( buf );
   }
   while( fread( &name_len, sizeof(int), 1, fp) ) {
      if( name_len > 200 )
         Rf_error( "One of the sequence names seems longer than 200 characters. "
           "Most likely this is not a valid BFA file." );
      fread( seq_name, sizeof(char), name_len, fp );
      fread( &seq_ori_len, sizeof(int), 1, fp );
      fread( &seq_len, sizeof(int), 1, fp );
      if( ( seq_ori_len >> 5 != seq_len ) && ( seq_ori_len >> 5 != seq_len - 1) )
         Rf_error( "Fields bfa.len and bfa_ori_len do not agree. This is not a "
           "valid BFA file." );
      fseek( fp, 2 * sizeof(bit64_t) * seq_len, SEEK_CUR);
      seqs.push_back( seq_meta_info( seq_ori_len, seq_name ) );
   }
   fclose( fp );
   
   SEXP res, names;
   PROTECT( res = allocVector( INTSXP, seqs.size() ) );   
   PROTECT( names = allocVector( STRSXP, seqs.size() ) );   
   int i = 0;
   for( std::deque< seq_meta_info >::iterator a = seqs.begin(); 
         a != seqs.end(); a++, i++ ) {
      INTEGER(res)[i] = a->len;      
      SET_STRING_ELT( names, i, mkChar( a->name.c_str() ) );
   }
   namesgets( res, names);   
   UNPROTECT(2);
   return res;
}   

