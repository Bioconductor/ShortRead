#include <stdlib.h>
#include "ShortRead.h"
#include "call.h"

/*
SIMU_0001_00000081/1	TGTACAGTATGTGAAGAGATTTGTTCTGAACCAAA	hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh	1	a	35	+	refseq	2210	0
*/

static const char *ELT_NMS[] = {
    "id", "sread", "quality", "nEquallyBestHits", "pairedEnd",
    "alignedLength", "strand", "chromosome", "position", "typeOfHit",
    "hitDetail"
};

static const int N_ELTS = sizeof(ELT_NMS) / sizeof(const char *);

int _read_soap(const char *fname, const char *csep, const char *commentChar,
               MARK_FIELD_FUNC * mark_func, SEXP ref, int offset)
{
    const int N_FIELDS = N_ELTS;
    gzFile file;
    char linebuf[LINEBUF_SIZE],
	**elt = (char **) R_alloc(N_FIELDS, sizeof(char*));
    int lineno = 0;

    file = _fopen(fname, "rb");

    _XSnap id = VECTOR_ELT(ref, 0),
        sread = VECTOR_ELT(ref, 1), quality = VECTOR_ELT(ref, 2);
    SEXP pairedEnd = VECTOR_ELT(ref, 4),
        chromosome = VECTOR_ELT(ref, 7), hitDetail = VECTOR_ELT(ref, 10);
    int *nEquallyBestHits = INTEGER(VECTOR_ELT(ref, 3)),
        *alignedLength = INTEGER(VECTOR_ELT(ref, 5)),
        *strand = INTEGER(VECTOR_ELT(ref, 6)),
        *position = INTEGER(VECTOR_ELT(ref, 8)),
        *typeOfHit = INTEGER(VECTOR_ELT(ref, 9));

    while (gzgets(file, linebuf, LINEBUF_SIZE) != NULL) {

        if (_linebuf_skip_p(linebuf, file, fname, lineno, commentChar)) {
            lineno++;
            continue;
        }

        /* field-ify */
        elt[0] = linebuf;
        for (int i = 1; i < N_FIELDS; ++i) {
            elt[i] = (*mark_func) (elt[i - 1], csep);
            if (elt[i] == elt[i - 1])
                error("too few fields, %s:%d", fname, lineno);
        }

        nEquallyBestHits[offset] = atoi(elt[3]);
        SET_STRING_ELT(pairedEnd, offset, mkChar(elt[4]));
        alignedLength[offset] = atoi(elt[5]);
        strand[offset] = _char_as_strand_int(*elt[6], fname, lineno);
        SET_STRING_ELT(chromosome, offset, mkChar(elt[7]));
        position[offset] = atoi(elt[8]);	/* leftmost-aligned, 1-based */
        typeOfHit[offset] = atoi(elt[9]);
        SET_STRING_ELT(hitDetail, offset, mkChar(elt[10]));
        /* 1-3: id, strand, quality */
        _APPEND_XSNAP(id, elt[0]);
        if (strand[offset] == 2) {
            _reverseComplement(elt[1]);
            _reverse(elt[2]);
        }
        _APPEND_XSNAP(sread, elt[1]);
        _APPEND_XSNAP(quality, elt[2]);

        lineno++;
        offset++;
    }
    return offset;
}

SEXP _AlignedRead_SOAP_make(SEXP ref, const char *qtype)
{
    SEXP s, t, nmspc = PROTECT(_get_namespace("ShortRead"));

    SEXP sfq;
    NEW_CALL(s, t, qtype, nmspc, 2);
    CSET_CDR(t, "quality", VECTOR_ELT(ref, 2));
    CEVAL_TO(s, nmspc, sfq);
    PROTECT(sfq);

    SEXP adf;
    NEW_CALL(s, t, ".SOAP_AlignedDataFrame", nmspc, 6);
    CSET_CDR(t, "nEquallyBestHits", VECTOR_ELT(ref, 3));
    CSET_CDR(t, "pairedEnd", VECTOR_ELT(ref, 4));
    CSET_CDR(t, "alignedLength", VECTOR_ELT(ref, 5));
    CSET_CDR(t, "typeOfHit", VECTOR_ELT(ref, 9));
    CSET_CDR(t, "hitDetail", VECTOR_ELT(ref, 10));
    CEVAL_TO(s, nmspc, adf);
    PROTECT(adf);

    SEXP aln;
    NEW_CALL(s, t, "AlignedRead", nmspc, 8);
    CSET_CDR(t, "sread", VECTOR_ELT(ref, 1));
    CSET_CDR(t, "id", VECTOR_ELT(ref, 0));
    CSET_CDR(t, "quality", sfq);
    CSET_CDR(t, "chromosome", VECTOR_ELT(ref, 7));
    CSET_CDR(t, "position", VECTOR_ELT(ref, 8));
    CSET_CDR(t, "strand", VECTOR_ELT(ref, 6));
    /* alignQuality */
    CSET_CDR(t, "alignData", adf);
    CEVAL_TO(s, nmspc, aln);

    UNPROTECT(3);
    return aln;
}

SEXP read_soap(SEXP files, SEXP qualityType, SEXP sep, SEXP commentChar)
{
    if (!IS_CHARACTER(files))
        Rf_error("'%s' must be '%s'", "files", "character()");
    if (!IS_CHARACTER(sep) || LENGTH(sep) != 1)
        Rf_error("'%s' must be '%s'", "sep", "character(1)");
    /* FIXME: !nzchar(sep[1]) */
    if (!IS_CHARACTER(commentChar) || LENGTH(commentChar) != 1)
        Rf_error("'%s' must be '%s'", "commentChar", "character(1)");
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
    SET_VECTOR_ELT(ref, 0, _NEW_XSNAP(nrec, "BString"));
    SET_VECTOR_ELT(ref, 1, _NEW_XSNAP(nrec, "DNAString"));
    SET_VECTOR_ELT(ref, 2, _NEW_XSNAP(nrec, "BString"));
    SET_VECTOR_ELT(ref, 3, NEW_INTEGER(nrec));	/* nEquallyBestHits */
    SET_VECTOR_ELT(ref, 4, NEW_STRING(nrec));	/* pairedEnd */
    SET_VECTOR_ELT(ref, 5, NEW_INTEGER(nrec));	/* alignedLength */
    SET_VECTOR_ELT(ref, 6, NEW_INTEGER(nrec));	/* strand */
    SET_VECTOR_ELT(ref, 7, NEW_STRING(nrec));	/* chromosome */
    SET_VECTOR_ELT(ref, 8, NEW_INTEGER(nrec));	/* position */
    SET_VECTOR_ELT(ref, 9, NEW_INTEGER(nrec));	/* typeOfHit */
    SET_VECTOR_ELT(ref, 10, NEW_STRING(nrec));	/* hitDetail */

    SEXP names = PROTECT(NEW_CHARACTER(N_ELTS));
    for (int i = 0; i < N_ELTS; ++i)
        SET_STRING_ELT(names, i, mkChar(ELT_NMS[i]));
    SET_ATTR(ref, R_NamesSymbol, names);
    UNPROTECT(1);

    const char *csep = translateChar(STRING_ELT(sep, 0));
    MARK_FIELD_FUNC *sep_func;  /* how to parse fields; minor efficiency */
    if (csep[0] != '\0' && csep[1] == '\0')
        sep_func = _mark_field_1;
    else
        sep_func = _mark_field_n;

    nrec = 0;
    for (int i = 0; i < LENGTH(files); ++i) {
        R_CheckUserInterrupt();
        nrec += _read_soap(CHAR(STRING_ELT(files, i)), csep,
                           CHAR(STRING_ELT(commentChar, 0)),
                           sep_func, ref, nrec);
    }
    _XSNAP_ELT(ref, 0);
    _XSNAP_ELT(ref, 1);
    _XSNAP_ELT(ref, 2);

    SEXP strand_lvls = PROTECT(_get_strand_levels());
    _as_factor_SEXP(VECTOR_ELT(ref, 6), strand_lvls);

    SEXP aln = _AlignedRead_SOAP_make(ref, qtype);
    UNPROTECT(2);
    return aln;
}
