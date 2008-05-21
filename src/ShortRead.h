#include <Rdefines.h>
#include "Biostrings_interface.h"

/* util.c */

typedef unsigned char (*DECODE_FUNC)(char); /* DNAdecode, RNAdecode */

DECODE_FUNC decoder(const char*);
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
SEXP alphabet_score(SEXP stringSet, SEXP vec);
SEXP alphabet_order(SEXP stringSet);
SEXP alphabet_duplicated(SEXP stringSet);
