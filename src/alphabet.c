#include "ShortRead.h"
#include <stdlib.h>

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
        Rf_error("'alphabet' must be character()");
    /* allocate and initialize the answer matrix */
    const int nrow = LENGTH(alphabet), ncol = INTEGER(width)[0];
    SEXP ans, dimnms, nms;
    PROTECT(ans = allocMatrix(INTSXP, nrow, ncol));
    PROTECT(dimnms = NEW_LIST(2));
    SET_VECTOR_ELT(dimnms, 0, alphabet);
    /* FIXME: Cycle dimnames? */
    PROTECT(nms = NEW_STRING(2));
    SET_STRING_ELT(nms, 0, mkChar("alphabet"));
    SET_STRING_ELT(nms, 1, mkChar("cycle"));
    SET_NAMES(dimnms, nms);
    SET_DIMNAMES(ans, dimnms);
    UNPROTECT(2);

    int *ansp = INTEGER(ans);   /* convenient pointer to data */
    memset(ansp, 0, LENGTH(ans) * sizeof(int)); /* initialize to 0 */

    /* set up a decoder for the string */
    const char *base = get_XStringSet_xsbaseclassname(stringSet);
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
     */
    cachedXStringSet cache = cache_XStringSet(stringSet);
    const int len = get_XStringSet_length(stringSet);
    for (i = 0; i < len; ++i) {
        cachedCharSeq seq = get_cachedXStringSet_elt(&cache, i);
        for (j = 0; j < seq.length; ++j) {
            int idx = map[decode(seq.seq[j])];
            if (idx >= 0)
                ansp[j * nrow + idx] += 1;
        }
    }

    UNPROTECT(1);
    return ans;
}

SEXP
alphabet_pair_by_cycle(SEXP stringSet1, SEXP stringSet2, SEXP width, SEXP alphabet1, SEXP alphabet2)
{
    const int MAX_MAP = 256;
    /* FIXME: check types of incoming arguments */
    if (get_XStringSet_length(stringSet1) != get_XStringSet_length(stringSet2))
        Rf_error("'stringSet1' and 'stringSet2' must have the same length");
    if (!IS_CHARACTER(alphabet1) || !IS_CHARACTER(alphabet2))
        Rf_error("'alphabet' must be list of character vectors");

    /* allocate and initialize the answer matrix */
    const int dim1 = LENGTH(alphabet1), dim2 = LENGTH(alphabet2), dim3 = INTEGER(width)[0];
    const int dim1xdim2 = dim1 * dim2;
    SEXP ans, dimnms, nms;
    PROTECT(ans = alloc3DArray(INTSXP, dim1, dim2, dim3));
    PROTECT(dimnms = NEW_LIST(3));
    SET_VECTOR_ELT(dimnms, 0, alphabet1);
    SET_VECTOR_ELT(dimnms, 1, alphabet2);
    /* FIXME: Cycle dimnames? */
    PROTECT(nms = NEW_STRING(3));
    SET_STRING_ELT(nms, 0, mkChar("base"));
    SET_STRING_ELT(nms, 1, mkChar("quality"));
    SET_STRING_ELT(nms, 3, mkChar("cycle"));
    SET_NAMES(dimnms, nms);
    SET_DIMNAMES(ans, dimnms);
    UNPROTECT(2);

    int *ansp = INTEGER(ans);   /* convenient pointer to data */
    memset(ansp, 0, LENGTH(ans) * sizeof(int)); /* initialize to 0 */

    /* set up a decoder for string1 and string2 */
    const char *base1 = get_XStringSet_xsbaseclassname(stringSet1);
    const char *base2 = get_XStringSet_xsbaseclassname(stringSet2);
    DECODE_FUNC decode1 = decoder(base1);
    DECODE_FUNC decode2 = decoder(base2);

    /* map between decoded character and offset into 'ans' */
    int i, j;
    int map1[MAX_MAP], map2[MAX_MAP];
    memset(map1, -1, MAX_MAP*sizeof(int)); /* default; ignore */
    memset(map2, -1, MAX_MAP*sizeof(int)); /* default; ignore */
    for (i = 0; i < LENGTH(alphabet1); ++i) {
        unsigned char c = (unsigned char) *CHAR(STRING_ELT(alphabet1, i));
        map1[c] = i;
    }
    for (i = 0; i < LENGTH(alphabet2); ++i) {
        unsigned char c = (unsigned char) *CHAR(STRING_ELT(alphabet2, i));
        map2[c] = i;
    }

    /* The main loop. Cache the string set for fast access, then
     * iterate over all strings, and over all characters in the
     * string. For each character, decode and map into the answer
     * matrix.
     *
     */
    cachedXStringSet cache1 = cache_XStringSet(stringSet1);
    cachedXStringSet cache2 = cache_XStringSet(stringSet2);
    const int len = get_XStringSet_length(stringSet1);
    for (i = 0; i < len; ++i) {
        cachedCharSeq seq1 = get_cachedXStringSet_elt(&cache1, i);
        cachedCharSeq seq2 = get_cachedXStringSet_elt(&cache2, i);
        for (j = 0; j < seq1.length; ++j) {
            int idx1 = map1[decode1(seq1.seq[j])];
            int idx2 = map2[decode2(seq2.seq[j])];
            if (idx1 >= 0 && idx2 >= 0)
                ansp[j * dim1xdim2 + idx2 * dim1 + idx1] += 1;
        }
    }

    UNPROTECT(1);
    return ans;
}

SEXP
alphabet_score(SEXP stringSet, SEXP score)
{
    /* FIXME: stringSet is XStringSet */
    const char *base = get_XStringSet_xsbaseclassname(stringSet);
    if (strcmp(base, "BString") != 0)
        Rf_error("'stringSet' must contain BString elements");
    if (!IS_NUMERIC(score) || LENGTH(score) != 256)
        Rf_error("'%s' must be '%s'", "score", "integer(256)");

    DECODE_FUNC decode = decoder(base);
    const int len = get_XStringSet_length(stringSet);
    int i, j;
    const double *dscore = REAL(score);

    SEXP ans;
    PROTECT(ans = NEW_NUMERIC(len));
    double *dans = REAL(ans);

    cachedXStringSet cache = cache_XStringSet(stringSet);
    for (i = 0; i < len; ++i) {
        cachedCharSeq seq = get_cachedXStringSet_elt(&cache, i);
        dans[i] = 0;
        for (j = 0; j < seq.length; ++j)
            dans[i] +=  dscore[decode(seq.seq[j])];
    }

    UNPROTECT(1);
    return ans;
}

SEXP
alphabet_as_int(SEXP stringSet, SEXP score)
{
    /* FIXME: stringSet is XStrinSet(1) or longer? */
    const char *base = get_XStringSet_xsbaseclassname(stringSet);
    if (strcmp(base, "BString") != 0)
        Rf_error("'stringSet' must contain BString elements");
    if (!IS_INTEGER(score) || LENGTH(score) != 256)
        Rf_error("'%s' must be '%s'", "score", "integer(256)");
    DECODE_FUNC decode = decoder(base);
    const int len = get_XStringSet_length(stringSet);

    cachedXStringSet cache = cache_XStringSet(stringSet);
    int i;

    cachedCharSeq seq = get_cachedXStringSet_elt(&cache, 0);
    int width = seq.length;
    int *ians;
    SEXP ans;
    for (i = 1; i < len && width > 0; ++i) {
        seq = get_cachedXStringSet_elt(&cache, i);
        if (seq.length != width) width = -1;
    }
    if (width >= 0) {           /* matrix */
        ans = PROTECT(allocMatrix(INTSXP, len, width));
        ians = INTEGER(ans);
    } else {                    /* list of int */
        ans = PROTECT(NEW_LIST(len));
    }

    const int *iscore = INTEGER(score);
    int j;
    for (i = 0; i < len; ++i) {
        seq = get_cachedXStringSet_elt(&cache, i);
        if (width >= 0) { /* int matrix */
            for (j = 0; j < seq.length; ++j)
                ians[len*j + i] =  iscore[decode(seq.seq[j])];
        } else {                /* list of ints */
            SET_VECTOR_ELT(ans, i, NEW_INTEGER(seq.length));
            ians = INTEGER(VECTOR_ELT(ans, i));
            for (j = 0; j < seq.length; ++j)
                ians[j] =  iscore[decode(seq.seq[j])];
        }
    }

    UNPROTECT(1);
    return ans;
}

/* rank / order / sort / duplicated */

typedef struct {
    int offset;
    cachedCharSeq ref;
} XSort;

int
compare_cachedCharSeq(const void *a, const void *b)
{
    const cachedCharSeq ra = ((const XSort*) a)->ref;
    const cachedCharSeq rb = ((const XSort*) b)->ref;

    const int diff = ra.length - rb.length;
    size_t len = diff < 0 ? ra.length : rb.length;
    int res = memcmp(ra.seq, rb.seq, len);
    return res == 0 ? diff : res;
}

void
_alphabet_order(cachedXStringSet cache, XSort *xptr, const int len)
{
    int i;

    for (i = 0; i < len; ++i) {
        xptr[i].offset=i;
        xptr[i].ref = get_cachedXStringSet_elt(&cache, i);
    }
    qsort(xptr, len, sizeof(XSort), compare_cachedCharSeq);
}

SEXP
alphabet_order(SEXP stringSet)
{
    /* FIXME: stringSet is XStringSet; non-zero len? */
    const int len = get_XStringSet_length(stringSet);
	if (len == 0)
		return NEW_INTEGER(0);
    cachedXStringSet cache = cache_XStringSet(stringSet);
    XSort *xptr = (XSort*) R_alloc(len, sizeof(XSort));
    _alphabet_order(cache, xptr, len);

    SEXP ans;
    PROTECT(ans = NEW_INTEGER(len));
    int *ians = INTEGER(ans);
    int i;
    for (i = 0; i < len; ++i)
        ians[i] = xptr[i].offset + 1;
    UNPROTECT(1);
    return ans;
}

SEXP
alphabet_duplicated(SEXP stringSet)
{
    /* FIXME: stringSet is XStringSet; non-zero len? */
    const int len = get_XStringSet_length(stringSet);
	if (len == 0)
		return NEW_LOGICAL(0);
    cachedXStringSet cache = cache_XStringSet(stringSet);
    XSort *xptr = (XSort*) R_alloc(len, sizeof(XSort));
    _alphabet_order(cache, xptr, len);

    SEXP ans;
    PROTECT(ans = NEW_LOGICAL(len));
    int *ians = INTEGER(ans);
    ians[xptr[0].offset]=0;
    int i;
    for (i = 1; i < len; ++i)
        ians[xptr[i].offset] = compare_cachedCharSeq(xptr+i-1, xptr+i) == 0;

    UNPROTECT(1);
    return ans;
}

SEXP
alphabet_rank(SEXP stringSet)
{
    /* integer vector of unique indices into sorted set */
    const int len = get_XStringSet_length(stringSet);
	if (len == 0)
		return NEW_INTEGER(0);
    cachedXStringSet cache = cache_XStringSet(stringSet);
    XSort *xptr = (XSort*) R_alloc(len, sizeof(XSort));
    _alphabet_order(cache, xptr, len);

    SEXP rank = PROTECT(NEW_INTEGER(len));
    int *irank = INTEGER(rank), i;
    irank[xptr[0].offset] = 1;
    for (i = 1; i < len; ++i) {
        if (compare_cachedCharSeq(&xptr[i-1], &xptr[i]) == 0) {
            irank[xptr[i].offset] = irank[xptr[i-1].offset];
        } else {
            irank[xptr[i].offset] = i + 1;
        }
    }

    UNPROTECT(1);
    return rank;
}

SEXP
aligned_read_rank(SEXP alignedRead, SEXP order, SEXP rho)
{
	if (LENGTH(order) == 0)
		return NEW_INTEGER(0);
    SEXP chr, str, pos, sread;
    PROTECT(chr = _get_SEXP(alignedRead, rho, "chromosome"));
    PROTECT(str = _get_SEXP(alignedRead, rho, "strand"));
    PROTECT(pos = _get_SEXP(alignedRead, rho, "position"));
    PROTECT(sread = _get_SEXP(alignedRead, rho, "sread"));
    int *c = INTEGER(chr), *s = INTEGER(str), *p = INTEGER(pos),
        *o = INTEGER(order), len = LENGTH(order);
    cachedXStringSet cache = cache_XStringSet(sread);
    XSort *xptr = (XSort*) R_alloc(2, sizeof(XSort));
    SEXP rank;
    PROTECT(rank = NEW_INTEGER(len));
    int *r = INTEGER(rank), i;
    xptr[0].ref = get_cachedXStringSet_elt(&cache, 0);
	r[o[0]-1] = 1;
	for (i = 1; i < len; ++i) {
		const int this = o[i]-1, prev=o[i-1]-1;
		xptr[i%2].ref = get_cachedXStringSet_elt(&cache, this);
		if (c[this] != c[prev] || s[this] != s[prev] ||
			p[this] != p[prev] || 
			compare_cachedCharSeq(xptr, xptr+1) != 0)
			r[this] = i + 1;
		else
			r[this] = r[prev];
	}
    UNPROTECT(5);
    return rank;
}
