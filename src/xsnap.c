/* 
 * An _XSNAP is a SEXP that contains sufficient information to create
 * (`snap') an XStringSet object from its content, without copying the
 * sequence data. It is allocated once to an initial size, and grows as
 * needed (the design goal is only one growth, using a simple
 * heuristic). Any `extra' allocation is not recovered at the end, but
 * carried forward until the DNAStringSet is garbage collected; the
 * design goal is that the extra allocation is no more than 10% of the
 * total object size.
 *
 * Basic usage is
 *
 *    SEXP lst = PROTECT(NEW_LIST(4));
 *    SET_VECTOR_ELT(lst, 0, _NEW_XSNAP(final_length_of_DNAStringSet));
 *    _APPEND_XSNAP(VECTOR_ELT(lst, 0), "ACTAGAC");
 *    SEXP xStringSet = PROTECT(_XSNAP_ELT(lst, 0, "DNAString"));
 *    UNPROTECT(2);
 *
 */

#include "ShortRead.h"

static const int INITIAL_SIZE = 32000;
static const double SCALE = 1.1;

_XSnap
_NEW_XSNAP(int nelt)
{
	_XSnap snap = PROTECT(NEW_LIST(5));

	SET_VECTOR_ELT(snap, 0, NEW_RAW(INITIAL_SIZE)); /* data */
	SET_VECTOR_ELT(snap, 1, NEW_INTEGER(nelt)); /* starts */
	SET_VECTOR_ELT(snap, 2, NEW_INTEGER(nelt)); /* widths */
	SET_VECTOR_ELT(snap, 3, NEW_INTEGER(2));	/* i, nelt */
	INTEGER(VECTOR_ELT(snap, 3))[0] = 0;
	INTEGER(VECTOR_ELT(snap, 3))[1] = nelt;
	SET_VECTOR_ELT(snap, 4, NEW_INTEGER(1));	/* offset */
	INTEGER(VECTOR_ELT(snap, 4))[0] = 0;
	UNPROTECT(1);
	return snap;
}

void
_APPEND_XSNAP(_XSnap snap, const char *str)
{
	const int width = strlen(str);

	const int i = INTEGER(VECTOR_ELT(snap, 3))[0];
	const int offset = INTEGER(VECTOR_ELT(snap, 4))[0];
	char *dest = (char *) RAW(VECTOR_ELT(snap, 0));

	int rlength = LENGTH(VECTOR_ELT(snap, 0));
	if (offset + width > rlength) {
		/* re-allocate */
		int nelt = INTEGER(VECTOR_ELT(snap, 3))[1];
		int sz = SCALE * rlength * nelt / i;
		SEXP buf = PROTECT(NEW_RAW(sz));
		memcpy((char *) RAW(buf), dest, offset);
		SET_VECTOR_ELT(snap, 0, buf);
		dest = (char *) RAW(buf);
		UNPROTECT(1);
	} 

	memcpy(dest + offset, str, width);
	INTEGER(VECTOR_ELT(snap, 1))[i] = offset + 1; /* R is 1-based */
	INTEGER(VECTOR_ELT(snap, 2))[i] = width;
	INTEGER(VECTOR_ELT(snap, 3))[0] = i + 1;
	INTEGER(VECTOR_ELT(snap, 4))[0] = offset + width;
}

SEXP
_XSnap_to_XStringSet(_XSnap snap, const char *baseclass)
{
	SEXP seq = VECTOR_ELT(snap, 0);
	int length = INTEGER(VECTOR_ELT(snap, 3))[0];
	int width = INTEGER(VECTOR_ELT(snap, 4))[0], i;

	ENCODE_FUNC encode = encoder(baseclass);
	char *str = (char *) RAW(seq);
	for (i = 0; i < width; ++i)
		str[i] = encode(str[i]);
	if (width < LENGTH(seq))	/* pad */
		memset(str + width, encode('N'), LENGTH(seq) - width);
	SEXP ptr = PROTECT(new_SharedVector("SharedRaw", seq));
	SEXP xstring = PROTECT(new_XVector(baseclass, ptr, 0, width));

	const int XSETCLASS_BUF = 40;
	if (strlen(baseclass) > XSETCLASS_BUF - 4)
		error("ShortRead internal error: *Set buffer too small");

	if (LENGTH(VECTOR_ELT(snap, 1)) != length) {
		SEXP tmp;
		tmp = VECTOR_ELT(snap, 1);
		SET_VECTOR_ELT(snap, 1, SET_LENGTH(tmp, length));
		tmp = VECTOR_ELT(snap, 2);
		SET_VECTOR_ELT(snap, 2, SET_LENGTH(tmp, length));
	}
	SEXP irange = 
		PROTECT(new_IRanges("IRanges", VECTOR_ELT(snap, 1),
							VECTOR_ELT(snap, 2), R_NilValue));
	SEXP xstringset = PROTECT(new_XStringSet(NULL, xstring, irange));
	UNPROTECT(4);
	return xstringset;
}

void
_XSNAP_ELT(SEXP x, int elt, const char *baseclass)
{
	SEXP xstringset = 
		PROTECT(_XSnap_to_XStringSet(VECTOR_ELT(x, elt), baseclass));
	SET_VECTOR_ELT(x, elt, xstringset);
	UNPROTECT(1);
}
