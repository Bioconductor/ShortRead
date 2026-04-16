#include <stdio.h>
#include <errno.h>
#include <deque>
#include <string>
#include <R.h>
#include <Rinternals.h>
#include "maqmap_m.h"

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
   if( (! Rf_isString(bfa_filename) ) || ( Rf_length(bfa_filename) != 1 ) )
      Rf_error( "First argument invalid: should be the filename." );

   fp = fopen( CHAR(STRING_ELT(bfa_filename,0)), "r" );
   if( !fp ) {
      char buf[300];
      snprintf( buf, 300, "Failed to open file '%s': %s (errno=%d)",
         CHAR(STRING_ELT(bfa_filename,0)), strerror(errno), errno );
      Rf_error( "%s", buf );
   }
/*   while( fread( &name_len, sizeof(int), 1, fp) ) {
      if( name_len > 200 )
         Rf_error( "sequence name >200 characters; invalid BFA file?" );
      (void) fread( seq_name, sizeof(char), name_len, fp );
      (void) fread( &seq_ori_len, sizeof(int), 1, fp );
      (void) fread( &seq_len, sizeof(int), 1, fp );
      if( ( seq_ori_len >> 5 != seq_len ) && ( seq_ori_len >> 5 != seq_len - 1) )
         Rf_error( "Fields bfa.len and bfa_ori_len do not agree. This is not a "
           "valid BFA file." );
      fseek( fp, 2 * sizeof(bit64_t) * seq_len, SEEK_CUR);
      seqs.push_back( seq_meta_info( seq_ori_len, seq_name ) );
   } */
    while (fread(&name_len, sizeof(int), 1, fp)) {
     if (name_len > 200)
        Rf_error("sequence name >200 characters; invalid BFA file?");
     if (fread(seq_name, sizeof(char), name_len, fp) != (size_t) name_len)
        Rf_error("truncated BFA file while reading sequence name");
     if (fread(&seq_ori_len, sizeof(int), 1, fp) != 1)
        Rf_error("truncated BFA file while reading seq_ori_len");
     if (fread(&seq_len, sizeof(int), 1, fp) != 1)
        Rf_error("truncated BFA file while reading seq_len");
     if ((seq_ori_len >> 5 != seq_len) && (seq_ori_len >> 5 != seq_len - 1))
        Rf_error("Fields bfa.len and bfa_ori_len do not agree. "
                 "This is not a valid BFA file.");
    fseek(fp, 2 * sizeof(bit64_t) * seq_len, SEEK_CUR);
    seqs.push_back(seq_meta_info(seq_ori_len, seq_name));
   }
   fclose( fp );
   
   SEXP res, names;
   PROTECT( res = Rf_allocVector( INTSXP, seqs.size() ) );   
   PROTECT( names = Rf_allocVector( STRSXP, seqs.size() ) );   
   int i = 0;
   for( std::deque< seq_meta_info >::iterator a = seqs.begin(); 
         a != seqs.end(); a++, i++ ) {
      INTEGER(res)[i] = a->len;      
      SET_STRING_ELT( names, i, Rf_mkChar( a->name.c_str() ) );
   }
   Rf_namesgets( res, names);   
   UNPROTECT(2);
   return res;
}   

