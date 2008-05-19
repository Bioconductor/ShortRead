#include "ShortRead.h"

/*
 * visit all sequences in a set, tallying character frequency as a
 * function of nucleotide position in the read.
 */
SEXP
alphabet_by_cycle(SEXP stringSet, SEXP width, SEXP alphabet)
{
    const int MAX_MAP = 256;
    /* FIXME: check types of incoming arguments */
    if (!IS_INTEGER(width) || LENGTH(width) != 1)
        Rf_error("'width' must be integer(1)");
    if (!IS_CHARACTER(alphabet))
        Rf_error("'alphabet' must be character");

    /* allocate and initialize the answer matrix */
    const int nrow = LENGTH(alphabet), ncol = INTEGER(width)[0];
    SEXP ans, nms;
    PROTECT(ans = allocMatrix(INTSXP, nrow, ncol));
    PROTECT(nms = NEW_LIST(2));
    SET_VECTOR_ELT(nms, 0, alphabet);
    setAttrib(ans, R_DimNamesSymbol, nms);
    UNPROTECT(1);

    int *ansp = INTEGER(ans);   /* convenient pointer to data */
    memset(ansp, 0, LENGTH(ans) * sizeof(int)); /* initialize to 0 */

    /* set up a decoder for the string */
    const char *base = get_XStringSet_baseClass(stringSet);
    DECODE_FUNC decode = decoder(base);

    /* map between decoded character and offset into 'ans' */
    int i, j;
    int map[MAX_MAP];
    memset(map, -1, MAX_MAP*sizeof(int)); /* default; ignore */
    for (i = 0; i < LENGTH(alphabet); ++i) {
        unsigned char c = (unsigned char) *CHAR(STRING_ELT(alphabet, i));
        map[c] = i;
    }    

    /* The main loop. Cache the string set for fast access, then
     * iterate over all strings, and over all characters in the
     * string. For each character, decode and map into the answer
     * matrix.
     *
     * FIXME: 
     */
    CachedXStringSet cache = new_CachedXStringSet(stringSet);
    const int len = get_XStringSet_length(stringSet);
    for (i = 0; i < len; ++i) {
        RoSeq seq = get_CachedXStringSet_elt_asRoSeq(&cache, i);
        for (j = 0; j < seq.nelt; ++j) {
            int idx = map[decode(seq.elts[j])];
            if (idx >= 0)
                ansp[j * nrow + idx] += 1;
        }
    }

    UNPROTECT(1);
    return ans;
}

SEXP
alphabet_score(SEXP stringSet, SEXP score)
{
    /* FIXME: stringSet is XStringSet */
    const char *base = get_XStringSet_baseClass(stringSet);
    if (strcmp(base, "BString") != 0)
	Rf_error("'stringSet' must contain BString elements");
    if (!IS_NUMERIC(score) || LENGTH(score) != 256)
	Rf_error("'score' must be numeric(256)");

    DECODE_FUNC decode = decoder(base);
    const int len = get_XStringSet_length(stringSet);
    int i, j;
    const double *dscore = REAL(score);

    SEXP ans;
    PROTECT(ans = NEW_NUMERIC(len));
    double *dans = REAL(ans);

    CachedXStringSet cache = new_CachedXStringSet(stringSet);
    for (i = 0; i < len; ++i) {
        RoSeq seq = get_CachedXStringSet_elt_asRoSeq(&cache, i);
	dans[i] = 0;
        for (j = 0; j < seq.nelt; ++j)
	    dans[i] +=  dscore[decode(seq.elts[j])];
    }

    UNPROTECT(1);
    return ans;
}
