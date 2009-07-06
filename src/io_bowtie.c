#include <stdlib.h>
#include "ShortRead.h"

/*
HWI-EAS88_1:1:1:83:277  -       chr1    163068612       AGAAGAATCCTTAAGGCTTGCTAGGCAGCAGTCTA     77777::::::::::::::::::::::::::::::     0       23
*/

static const char *ELT_NMS[] = {
    "id", "strand", "chromosome", "position", "sread", "quality",
    "similar", "mismatch"
};
static const int N_ELTS = sizeof(ELT_NMS) / sizeof(const char*);

int
_read_bowtie(const char *fname, const char *commentChar,
             SEXP ref, int offset)
{
    const int N_FIELDS = 8;
    gzFile *file;
    char linebuf[LINEBUF_SIZE], *elt[N_FIELDS];
    int lineno = 0;

    file = _fopen(fname, "rb");

	_XSnap id = VECTOR_ELT(ref, 0),
		sread = VECTOR_ELT(ref, 4),
		quality = VECTOR_ELT(ref, 5);
    SEXP chromosome = VECTOR_ELT(ref, 2),
        mismatch = VECTOR_ELT(ref, 7);
    int *strand = INTEGER(VECTOR_ELT(ref, 1)),
      *position = INTEGER(VECTOR_ELT(ref, 3)),
      *similar = INTEGER(VECTOR_ELT(ref, 6));

    while (gzgets(file, linebuf, LINEBUF_SIZE) != NULL) {

		if (*linebuf == *commentChar) {
            lineno++;
            continue;
        }
        lineno++;

		int n_fields = _mark_field_0(linebuf, elt, N_FIELDS);
		if (n_fields != N_FIELDS) {
			gzclose(file);
			error("incorrect number of fields (%d) %s:%d",
				  n_fields, fname, lineno);
		}

        _APPEND_XSNAP(id, elt[0]);
        strand[offset] = _char_as_strand_int(*elt[1], fname, lineno);
        SET_STRING_ELT(chromosome, offset, mkChar(elt[2]));
        position[offset] = atoi(elt[3]) + 1; /* leftmost-aligned, 0-based */
        if (strand[offset] == 1) {
            _reverseComplement(elt[4]);
            _reverse(elt[5]);
        }
		_APPEND_XSNAP(sread, elt[4]);
		_APPEND_XSNAP(quality, elt[5]);
        similar[offset] = atoi(elt[6]); /* previous: 'reserved' */
        SET_STRING_ELT(mismatch, offset, mkChar(elt[7]));
        offset++;
    }
    return offset;
}

#define NEW_CALL(S, T, NAME, ENV, N) \
    PROTECT(S = T = allocList(N)); \
    SET_TYPEOF(T, LANGSXP); \
    SETCAR(T, findFun(install(NAME), ENV)); \
    T = CDR(T)
#define CSET_CDR(T, NAME, VALUE) \
    SETCAR(T, VALUE); \
    SET_TAG(T, install(NAME)); \
    T = CDR(T)
#define CEVAL_TO(S, ENV, GETS) \
    GETS = eval(S, ENV); \
    UNPROTECT(1)

SEXP
_AlignedRead_Bowtie_make(SEXP ref, const char *qtype)
{
    SEXP s, t, nmspc = PROTECT(_get_namespace("ShortRead"));

    SEXP sfq;                   /* SFastqQuality by default */
    NEW_CALL(s, t, qtype, nmspc, 2);
    CSET_CDR(t, "quality", VECTOR_ELT(ref, 5));
    CEVAL_TO(s, nmspc, sfq);
    PROTECT(sfq);

    SEXP adf;
    NEW_CALL(s, t, ".Bowtie_AlignedDataFrame", nmspc, 3);
    CSET_CDR(t, "similar", VECTOR_ELT(ref, 6));
    CSET_CDR(t, "mismatch", VECTOR_ELT(ref, 7));
    CEVAL_TO(s, nmspc, adf);
    PROTECT(adf);

    SEXP aln;
    NEW_CALL(s, t, "AlignedRead", nmspc, 8);
    CSET_CDR(t, "id", VECTOR_ELT(ref, 0));
    CSET_CDR(t, "sread", VECTOR_ELT(ref, 4));
    CSET_CDR(t, "quality", sfq);
    CSET_CDR(t, "chromosome", VECTOR_ELT(ref, 2));
    CSET_CDR(t, "position", VECTOR_ELT(ref, 3));
    CSET_CDR(t, "strand", VECTOR_ELT(ref, 1));
    /* alignQuality */
    CSET_CDR(t, "alignData", adf);
    CEVAL_TO(s, nmspc, aln);

    UNPROTECT(3);
    return aln;
}

#undef NEW_CALL
#undef CSET_CDR
#undef CEVAL_TO

SEXP
read_bowtie(SEXP files, SEXP qualityType, SEXP sep, SEXP commentChar)
{
    if (!IS_CHARACTER(files))
        Rf_error("'files' must be 'character()'");
    if (!IS_CHARACTER(sep) || LENGTH(sep) != 1 || 
		*CHAR(STRING_ELT(sep, 0)) != '\t')
        Rf_error("'sep' must be '\t'"); 
    if (!IS_CHARACTER(commentChar) || LENGTH(commentChar) != 1)
        Rf_error("'commentChar' must be character(1)");
    if (LENGTH(STRING_ELT(commentChar, 0)) != 1)
        Rf_error("'nchar(commentChar[[1]])' must be 1 but is %d",
                 LENGTH(STRING_ELT(commentChar, 0)));
    if (!IS_CHARACTER(qualityType) || LENGTH(qualityType) != 1)
        Rf_error("'%s' must be '%s'", "qualityType", "character(1)");
    const char *qtype = CHAR(STRING_ELT(qualityType, 0));
    if (strcmp(qtype, "SFastqQuality") != 0 &&
        strcmp(qtype, "FastqQuality") != 0)
        Rf_error("'%s' must be '%s'", "qualityType",
                 "SFastqQuality' or 'FastqQuality");

    int nrec = _count_lines_sum(files);
    SEXP ref = PROTECT(NEW_LIST(N_ELTS));
	SET_VECTOR_ELT(ref, 0, _NEW_XSNAP(nrec));  /* id */
    SET_VECTOR_ELT(ref, 1, NEW_INTEGER(nrec)); /* strand */
    SET_VECTOR_ELT(ref, 2, NEW_STRING(nrec)); /* chromosome */
    SET_VECTOR_ELT(ref, 3, NEW_INTEGER(nrec)); /* position */
	SET_VECTOR_ELT(ref, 4, _NEW_XSNAP(nrec)); /* sread */
	SET_VECTOR_ELT(ref, 5, _NEW_XSNAP(nrec)); /* quality */
    SET_VECTOR_ELT(ref, 6, NEW_INTEGER(nrec)); /* similar */
    SET_VECTOR_ELT(ref, 7, NEW_STRING(nrec)); /* mismatch encoding */

    SEXP names = PROTECT(NEW_CHARACTER(N_ELTS));
    for (int i = 0; i < N_ELTS; ++i)
        SET_STRING_ELT(names, i, mkChar(ELT_NMS[i]));
    SET_ATTR(ref, R_NamesSymbol, names);
    UNPROTECT(1);

    nrec = 0;
    for (int i = 0; i < LENGTH(files); ++i) {
        R_CheckUserInterrupt();
        nrec += _read_bowtie(
            CHAR(STRING_ELT(files, i)),
            CHAR(STRING_ELT(commentChar, 0)),
            ref, nrec);
    }
    _XSNAP_ELT(ref, 0, "BString");
    _XSNAP_ELT(ref, 4, "DNAString");
    _XSNAP_ELT(ref, 5, "BString");

    SEXP strand_lvls = PROTECT(_get_strand_levels());
    _as_factor_SEXP(VECTOR_ELT(ref, 1), strand_lvls);
	UNPROTECT(1);

    SEXP aln = _AlignedRead_Bowtie_make(ref, qtype);

    UNPROTECT(1);
    return aln;
}
