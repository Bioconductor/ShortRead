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

/* _Snap */

struct _Snap {
	int i_entries, n_entries;
	char **entries;
};

struct _Snap *
snap_new(int n_entries)
{
	struct _Snap *snap = Calloc(1, struct _Snap);
	snap->entries = Calloc(n_entries, char *);
	for (int i = 0; i < n_entries; ++i)
		snap->entries[i] = NULL;
	snap->n_entries = n_entries;
	snap->i_entries = 0;
	return snap;
}

void
snap_encode(struct _Snap *snap, ENCODE_FUNC encode)
{
	for (int i = 0; i < snap->i_entries; ++i)
		for (char *c = snap->entries[i]; *c != '\0'; ++c)
			*c = encode(*c);
}

void
snap_free(struct _Snap *snap)
{
	Free(snap->entries);
	Free(snap);
}

/* _XSnap: wrap in ExternalPtr as interface and to ensure garbage
 * collection */

void
_xsnap_finalizer(SEXP xsnap)
{
	struct _Snap *snap = R_ExternalPtrAddr(xsnap);
	if (!snap) return;
	snap_free(snap);
	R_ClearExternalPtr(xsnap);
}

_XSnap
_NEW_XSNAP(int nelt)
{
	struct _Snap *snap = snap_new(nelt);
	SEXP xsnap = PROTECT(R_MakeExternalPtr(snap, mkString("XSnap"),
										   R_NilValue));
	R_RegisterCFinalizerEx(xsnap, _xsnap_finalizer, TRUE);
	UNPROTECT(1);
	return xsnap;
}

void
_APPEND_XSNAP(_XSnap xsnap, const char *str)
{
	struct _Snap *snap = (struct _Snap *) R_ExternalPtrAddr(xsnap);
	if (snap->i_entries >= snap->n_entries)
		Rf_error("ShortRead internal: too many 'snap' entries");
	snap->entries[snap->i_entries] = Calloc(strlen(str) + 1, char);
	strcpy(snap->entries[snap->i_entries], str);
	snap->i_entries += 1;
}

SEXP
_XSnap_to_XStringSet(_XSnap xsnap, const char *baseclass)
{
	struct _Snap *snap = (struct _Snap *) R_ExternalPtrAddr(xsnap);

	snap_encode(snap, encoder(baseclass));
	SEXP width = PROTECT(NEW_INTEGER(snap->i_entries));
	int i, *w = INTEGER(width);
	for (i = 0; i < snap->i_entries; ++i)
		w[i] = strlen(snap->entries[i]);
	
	const int XSETCLASS_BUF = 40;
	if (strlen(baseclass) > XSETCLASS_BUF - 4)
		Rf_error("ShortRead internal: *Set buffer too small");
	char class[XSETCLASS_BUF];
	sprintf(class, "%sSet", baseclass);

	SEXP ans = PROTECT(alloc_XRawList(class, baseclass, width));
	cachedXVectorList cached_ans = cache_XVectorList(ans);
	cachedCharSeq elt;

	for (i = 0; i < snap->i_entries; ++i) {
		elt = get_cachedXRawList_elt(&cached_ans, i);
		memcpy((char *) elt.seq, snap->entries[i], w[i]);
	}

	UNPROTECT(2);
	return ans;
}

void
_XSNAP_ELT(SEXP x, int elt, const char *baseclass)
{
	SEXP xstringset = 
		PROTECT(_XSnap_to_XStringSet(VECTOR_ELT(x, elt), baseclass));
	SET_VECTOR_ELT(x, elt, xstringset);
	UNPROTECT(1);
}
