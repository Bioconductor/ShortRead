#include <Rdefines.h>
#include "Biostrings_interface.h"

/* util.c */

typedef unsigned char (*DECODE_FUNC)(char); /* DNAdecode, RNAdecode */
DECODE_FUNC decoder(const char*);

typedef char * (MARK_FIELD_FUNC)(char *, const char *);
MARK_FIELD_FUNC _mark_field_1;	/* nchar(delim) == 1 */
MARK_FIELD_FUNC _mark_field_n;	/* nchar(delim) != 1 */

int _rtrim(char *linebuf);
void _solexa_to_IUPAC(char *linebuf);
SEXP count_lines(SEXP files);

/* io.c */

SEXP read_prb_as_character(SEXP file, SEXP cycles);
SEXP read_solexa_fastq(SEXP files);
SEXP read_XStringSet_columns(SEXP file, SEXP colIndex,
                             SEXP colClasses, SEXP delim,
                             SEXP header, SEXP commentChar);

/* alphabet.c */

SEXP alphabet_by_cycle(SEXP stringSet, SEXP width, SEXP alphabet);
SEXP alphabet_pair_by_cycle(SEXP stringSet1, SEXP stringSet2, SEXP width,
                            SEXP alphabet1, SEXP alphabet2);
SEXP alphabet_score(SEXP stringSet, SEXP score);
SEXP alphabet_as_int(SEXP stringSet, SEXP score);
SEXP alphabet_order(SEXP stringSet);
SEXP alphabet_duplicated(SEXP stringSet);
SEXP alphabet_rank(SEXP stringSet);

/* read_maq_map.c */

SEXP read_maq_map(SEXP filename, SEXP maxreads);

/* pileup.c */

SEXP pileup(SEXP start, SEXP fraglength, SEXP chrlength, SEXP dir,
            SEXP readlength, SEXP offset);

/* readBfaToc.c */

SEXP readBfaToc(SEXP bfa_filename);
