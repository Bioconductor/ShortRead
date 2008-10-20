#ifndef _SHORTREAD_H_
#define _SHORTREAD_H_

#ifdef __cplusplus
extern "C" {
#endif


#include <Rdefines.h>
#include "IRanges_interface.h"
#include "Biostrings_interface.h"

/* util.c */

typedef unsigned char (*DECODE_FUNC)(char); /* DNAdecode, RNAdecode */
DECODE_FUNC decoder(const char*);

SEXP _get_namespace(const char *pkg);
SEXP _get_strand_levels();

typedef char * (MARK_FIELD_FUNC)(char *, const char *);
MARK_FIELD_FUNC _mark_field_1;	/* nchar(delim) == 1 */
MARK_FIELD_FUNC _mark_field_n;	/* nchar(delim) != 1 */

int _rtrim(char *linebuf);
void _solexa_to_IUPAC(char *linebuf);
SEXP _CharAEAE_to_XStringSet(CharAEAE* aeae, const char *clsName);
void _as_factor_SEXP(SEXP vec, SEXP lvls);
void _as_factor(SEXP vec, const char **levels, const int n_lvls);
int _count_lines_sum(SEXP files);
SEXP count_lines(SEXP files);
SEXP _get_SEXP(SEXP from, SEXP rho, const char *with);

/* io.c */

SEXP read_prb_as_character(SEXP file, SEXP cycles);
SEXP read_solexa_fastq(SEXP files);
SEXP read_XStringSet_columns(SEXP file, SEXP colIndex,
                             SEXP colClasses, SEXP delim,
                             SEXP header, SEXP commentChar);
SEXP read_solexa_export(SEXP files, SEXP sep, SEXP commentChar);

/* alphabet.c */

SEXP alphabet_by_cycle(SEXP stringSet, SEXP width, SEXP alphabet);
SEXP alphabet_pair_by_cycle(SEXP stringSet1, SEXP stringSet2, SEXP width,
                            SEXP alphabet1, SEXP alphabet2);
SEXP alphabet_score(SEXP stringSet, SEXP score);
SEXP alphabet_as_int(SEXP stringSet, SEXP score);
SEXP alphabet_order(SEXP stringSet);
SEXP alphabet_duplicated(SEXP stringSet);
SEXP alphabet_rank(SEXP stringSet);
SEXP aligned_read_rank(SEXP stringSet, SEXP order, SEXP rho);

/* read_maq_map.c */

SEXP read_maq_map(SEXP filename, SEXP maxreads, SEXP maq_longread);

/* pileup.c */

SEXP pileup(SEXP start, SEXP fraglength, SEXP chrlength, SEXP dir,
            SEXP readlength, SEXP offset);

/* readBfaToc.c */

SEXP readBfaToc(SEXP bfa_filename);

#ifdef __cplusplus
}
#endif

#endif /* _SHORTREAD_H_ */
