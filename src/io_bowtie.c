#include <stdlib.h>
#include "ShortRead.h"

/*
HWI-EAS88_1:1:1:83:277  -       chr1    163068612       AGAAGAATCCTTAAGGCTTGCTAGGCAGCAGTCTA     77777::::::::::::::::::::::::::::::     0       23
*/

static const char *ELT_NMS[] = {
    "id", "strand", "chromosome", "position", "sread", "quality",
    "mismatch"
};
static const int N_ELTS = sizeof(ELT_NMS) / sizeof(const char*);

int
_read_bowtie(const char *fname, const char *csep,
             const char *commentChar,
             MARK_FIELD_FUNC *mark_func,
             SEXP ref, int offset,
             CharAEAE *sread, CharAEAE *quality)
{
    const int N_FIELDS = 8;
    gzFile *file;
    char linebuf[LINEBUF_SIZE], *elt[N_FIELDS];
    int lineno = 0;

    file = _fopen(fname, "rb");

    SEXP chromosome = VECTOR_ELT(ref, 2),
        mismatch = VECTOR_ELT(ref, 6);
    int *strand = INTEGER(VECTOR_ELT(ref, 1)),
        *position = INTEGER(VECTOR_ELT(ref, 3));

    while (gzgets(file, linebuf, LINEBUF_SIZE) != NULL) {

        if (_linebuf_skip_p(linebuf, file, fname, lineno,
                            commentChar)) {
            lineno++;
            continue;
        }

        /* field-ify */
        elt[0] = linebuf;
        for (int i = 1; i < N_FIELDS; ++i) {
            elt[i] = (*mark_func)(elt[i-1], csep);
            if (elt[i] == elt[i-1])
                error("too few fields, %s:%d", fname, lineno);
        }
        /* elt[0] (id) ignored */
        if (*elt[1] == '\0')
            strand[offset] = NA_INTEGER;
        else {
            switch (*elt[1]) {
            case '-':
                strand[offset] = 1;
                break;
            case '+':
                strand[offset] = 2;
                break;
            default:
                error("invalid 'strand' field '%s', %s:%d",
                      *elt[1], fname, lineno);
                break;
            }
        }
        SET_STRING_ELT(chromosome, offset, mkChar(elt[2]));
        position[offset] = atoi(elt[3]) + 1; /* input is leftmost-aligned, 0-based */
        if (strand[offset] == 1) {
            _reverseComplement(elt[4]);
            _reverse(elt[5]);
        }
        append_string_to_CharAEAE(sread, elt[4]);
        append_string_to_CharAEAE(quality, elt[5]);
        /* 'internal', ignored */
        SET_STRING_ELT(mismatch, offset, mkChar(elt[7]));
        lineno++;
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
    NEW_CALL(s, t, ".Bowtie_AlignedDataFrame", nmspc, 2);
    CSET_CDR(t, "mismatch", VECTOR_ELT(ref, 6));
    CEVAL_TO(s, nmspc, adf);
    PROTECT(adf);

    SEXP aln;
    SEXP strand_lvls = PROTECT(_get_strand_levels());
    _as_factor_SEXP(VECTOR_ELT(ref, 1), strand_lvls);
    NEW_CALL(s, t, "AlignedRead", nmspc, 7);
    CSET_CDR(t, "sread", VECTOR_ELT(ref, 4));
    CSET_CDR(t, "quality", sfq);
    CSET_CDR(t, "chromosome", VECTOR_ELT(ref, 2));
    CSET_CDR(t, "position", VECTOR_ELT(ref, 3));
    CSET_CDR(t, "strand", VECTOR_ELT(ref, 1));
    /* alignQuality */
    CSET_CDR(t, "alignData", adf);
    CEVAL_TO(s, nmspc, aln);

    UNPROTECT(4);
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
    if (!IS_CHARACTER(sep) || LENGTH(sep) != 1)
        Rf_error("'sep' must be character(1)"); 
    /* FIXME: !nzchar(sep[1]) */
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
    /* SET_VECTOR_ELT(ref, 0, NEW_STRING(nrec)); id, ignored */
    SET_VECTOR_ELT(ref, 1, NEW_INTEGER(nrec)); /* strand */
    SET_VECTOR_ELT(ref, 2, NEW_STRING(nrec)); /* chromosome */
    SET_VECTOR_ELT(ref, 3, NEW_INTEGER(nrec)); /* position */
    CharAEAE
        sread = new_CharAEAE(nrec, 0),
        quality = new_CharAEAE(nrec, 0);
    /* 'reserved'; ignore */
    SET_VECTOR_ELT(ref, 6, NEW_STRING(nrec)); /* mismatch encoding */

    SEXP names = PROTECT(NEW_CHARACTER(N_ELTS));
    for (int i = 0; i < N_ELTS; ++i)
        SET_STRING_ELT(names, i, mkChar(ELT_NMS[i]));
    SET_ATTR(ref, R_NamesSymbol, names);
    UNPROTECT(1);

    const char *csep = translateChar(STRING_ELT(sep, 0));
    MARK_FIELD_FUNC *sep_func;/* how to parse fields; minor efficiency */
    if (csep[0] != '\0' && csep[1] == '\0')
        sep_func = _mark_field_1;
    else
        sep_func = _mark_field_n;

    nrec = 0;
    for (int i = 0; i < LENGTH(files); ++i) {
        R_CheckUserInterrupt();
        nrec += _read_bowtie(
            CHAR(STRING_ELT(files, i)), csep,
            CHAR(STRING_ELT(commentChar, 0)),
            sep_func, ref, nrec, &sread, &quality);
    }
    SET_VECTOR_ELT(ref, 4, _CharAEAE_to_XStringSet(&sread, "DNAString"));
    SET_VECTOR_ELT(ref, 5, _CharAEAE_to_XStringSet(&quality, "BString"));
    SEXP aln = _AlignedRead_Bowtie_make(ref, qtype);
    UNPROTECT(1);
    return aln;
}
