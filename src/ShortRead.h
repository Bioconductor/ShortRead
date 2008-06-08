#include <Rdefines.h>
#include "Biostrings_interface.h"

/* util.c */

typedef unsigned char (*DECODE_FUNC)(char); /* DNAdecode, RNAdecode */

DECODE_FUNC decoder(const char*);
char *_mark_field(char *curr, const char* delim);
int _rtrim(char *linebuf);
void _solexa_to_IUPAC(char *linebuf);
SEXP count_lines(SEXP files);

/* io.c */

SEXP read_solexa_fastq(SEXP files);
SEXP read_XStringSet_columns(SEXP file, SEXP colIndex,
                             SEXP colClasses, SEXP delim,
                             SEXP header);

/* alphabet.c */

SEXP alphabet_by_cycle(SEXP stringSet, SEXP width, SEXP alphabet);
SEXP alphabet_score(SEXP stringSet, SEXP score);
SEXP alphabet_as_int(SEXP stringSet, SEXP score);
SEXP alphabet_order(SEXP stringSet);
SEXP alphabet_duplicated(SEXP stringSet);
SEXP alphabet_rank(SEXP stringSet);

/* read_maq_map.c */

SEXP read_maq_map( SEXP filename, SEXP maxreads );
