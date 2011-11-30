#ifndef _SHORTREAD_H_
#define _SHORTREAD_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <zlib.h>
#include <Rdefines.h>
#include "IRanges_interface.h"
#include "Biostrings_interface.h"

/* util.c */

    typedef unsigned char (*DECODE_FUNC) (char);	/* DNAdecode, RNAdecode */
    typedef char (*ENCODE_FUNC) (char);	/* DNAdecode, RNAdecode */
    DECODE_FUNC decoder(const char *);
    ENCODE_FUNC encoder(const char *);

    void _reverse(char *);
    void _reverseComplement(char *);

    SEXP _get_namespace(const char *pkg);
    SEXP _get_strand_levels();
    int _char_as_strand_int(const char c, const char *fname, const int lineno);

    typedef char *(MARK_FIELD_FUNC) (char *, const char *);
    MARK_FIELD_FUNC _mark_field_1;	/* nchar(delim) == 1 */
    MARK_FIELD_FUNC _mark_field_n;	/* nchar(delim) != 1 */
    int _mark_field_0(char *, char **, const int);

    extern const int LINEBUF_SIZE;
    gzFile *_fopen(const char *, const char *);
    int _linebuf_skip_p(char *, gzFile *, const char *, int, const char *);

    int _rtrim(char *linebuf);
    void _solexa_to_IUPAC(char *linebuf);
    void _as_factor_SEXP(SEXP vec, SEXP lvls);
    void _as_factor(SEXP vec, const char **levels, const int n_lvls);
    int _count_lines_sum(SEXP files);
    SEXP count_lines(SEXP files);
    SEXP count_ipar_int_recs(SEXP files);
    SEXP _get_SEXP(SEXP from, SEXP rho, const char *with);

/* xstring_util.c */

    typedef SEXP _XSnap;
    _XSnap _NEW_XSNAP(int nelt, const char *baseclass);
    void _APPEND_XSNAP(_XSnap snap, const char *str);
    void _XSNAP_ELT(SEXP x, int elt);

/* io.c */

    SEXP write_fastq(SEXP id, SEXP sread, SEXP quality,
                     SEXP fname, SEXP fmode, SEXP full, SEXP max_width);
    SEXP read_prb_as_character(SEXP file, SEXP asSolexa);
    SEXP read_solexa_fastq(SEXP files, SEXP withIds);
    SEXP read_XStringSet_columns(SEXP files, SEXP header, SEXP sep,
                                 SEXP colIndex, SEXP colClasses,
                                 SEXP nrows, SEXP skip, SEXP commentChar);
    SEXP read_solexa_export(SEXP files, SEXP sep, SEXP commentChar,
                            SEXP withFlags);

/* io_bowtie.c, io_soap.c */

    SEXP read_bowtie(SEXP files, SEXP qualityType, SEXP sep, SEXP commentChar);
    SEXP read_soap(SEXP files, SEXP qualityType, SEXP sep, SEXP commentChar);

/* alphabet.c */

    SEXP alphabet_by_cycle(SEXP stringSet, SEXP width, SEXP alphabet);
    SEXP alphabet_pair_by_cycle(SEXP stringSet1, SEXP stringSet2, SEXP width,
                                SEXP alphabet1, SEXP alphabet2);
    SEXP alphabet_score(SEXP stringSet, SEXP score);
    SEXP alphabet_as_int(SEXP stringSet, SEXP score);
    SEXP alphabet_order(SEXP stringSet);
    SEXP alphabet_duplicated(SEXP stringSet);
    SEXP alphabet_rank(SEXP stringSet);
    SEXP aligned_read_rank(SEXP stringSet, SEXP order, SEXP withSread,
                           SEXP rho);

/* read_maq_map.c */

    SEXP read_maq_map(SEXP filename, SEXP maxreads, SEXP maq_longread);

/* pileup.c */

    SEXP pileup(SEXP start, SEXP fraglength, SEXP chrlength, SEXP dir,
                SEXP readlength, SEXP offset);

/* readBfaToc.c */

    SEXP readBfaToc(SEXP bfa_filename);

/* sampler */

    SEXP sampler_new(SEXP n);
    SEXP sampler_add(SEXP s, SEXP bin);
    SEXP sampler_status(SEXP s);
    SEXP sampler_as_XStringSet(SEXP s);

    SEXP streamer_new(SEXP n);
    SEXP streamer_add(SEXP s, SEXP bin);
    SEXP streamer_status(SEXP s);
    SEXP streamer_as_XStringSet(SEXP s);

#ifdef __cplusplus
}
#endif
#endif                          /* _SHORTREAD_H_ */
