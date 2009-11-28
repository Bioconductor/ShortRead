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
 *    SET_VECTOR_ELT(lst, 0, _NEW_XSNAP(final_length, "DNAString"));
 *    _APPEND_XSNAP(VECTOR_ELT(lst, 0), "ACTAGAC");
 *    SEXP xStringSet = PROTECT(_XSNAP_ELT(lst, 0));
 *    UNPROTECT(2);
 *
 */

#include "ShortRead.h"

/* 50m 2x100 reads; 10G reads */
static const int GROUP_SIZE = 100000000;
static const int N_GROUP = 1000;

enum { GROUPS_IDX = 0, READS_PER_GROUP_IDX, GROUP_IDX, START_IDX, 
	   WIDTH_IDX, STATE_IDX, BASENAME_IDX };

#define GET_GRP_I(snap, i) \
	VECTOR_ELT(VECTOR_ELT((snap), GROUPS_IDX), i)
#define SET_GRP_I(snap, i, value) \
	SET_VECTOR_ELT(VECTOR_ELT((snap), GROUPS_IDX), i, (value))
#define GET_CURR_GRP(snap) GET_GRP_I((snap), GET_STATE((snap))[1])
#define SET_CURR_GRP(snap, value) \
	SET_GRP_I((snap), GET_STATE((snap))[1], (value))
#define GET_STATE(snap) INTEGER(VECTOR_ELT((snap), STATE_IDX))

void
_trim_group_i(_XSnap snap, int i, int len)
{
	SEXP grp = GET_GRP_I(snap, i);
	SET_LENGTH(grp, len);
	SET_GRP_I(snap, i, grp);
}

_XSnap
_NEW_XSNAP(int nelt, const char *baseclass)
{
	_XSnap snap = PROTECT(NEW_LIST(7));

	SET_VECTOR_ELT(snap, GROUPS_IDX, NEW_LIST(N_GROUP)); /* GROUP list */
	SET_VECTOR_ELT(snap, READS_PER_GROUP_IDX, 
				   NEW_INTEGER(GROUP_SIZE)); /* group sizes */
	SET_VECTOR_ELT(snap, GROUP_IDX, NEW_INTEGER(nelt)); /* group */
	SET_VECTOR_ELT(snap, START_IDX, NEW_INTEGER(nelt)); /* starts */
	SET_VECTOR_ELT(snap, WIDTH_IDX, NEW_INTEGER(nelt)); /* widths */
	SET_VECTOR_ELT(snap, STATE_IDX, NEW_INTEGER(5));	/* state */
	SET_VECTOR_ELT(snap, BASENAME_IDX, NEW_CHARACTER(1));
	int *state = GET_STATE(snap);
	state[0] = 0;							/* i */
	state[1] = 0;							/* current group */
	state[2] = 0;							/* offset */
	state[3] = nelt;						/* nelt */
	state[4] = N_GROUP;						/* ngroups */
	SET_CURR_GRP(snap, NEW_RAW(GROUP_SIZE));
	SET_STRING_ELT(VECTOR_ELT(snap, BASENAME_IDX), 0, mkChar(baseclass));
	UNPROTECT(1);
	return snap;
}

void
_APPEND_XSNAP(_XSnap snap, const char *str)
{
	const int width = strlen(str);

	int *state = GET_STATE(snap);
	SEXP grp = GET_CURR_GRP(snap);
	if (state[2] + width > LENGTH(grp)) { /* grow */
		if (state[1] >= state[4]) {	  /* grow groups */
			int new_sz = state[4] * 2;
			SEXP grps = NEW_LIST(state[4]); /* no PROTECT needed here */
			for (int i = 0; i < state[1]; ++i)
				SET_VECTOR_ELT(grps, i, GET_GRP_I(snap, i));
			for (int i = state[i]; i < new_sz; ++i)
				SET_VECTOR_ELT(grps, i, R_NilValue);
			SET_VECTOR_ELT(snap, GROUPS_IDX, grps);
			int *orig_reads_per_group = 
				INTEGER(VECTOR_ELT(snap, READS_PER_GROUP_IDX));
			SEXP reads_per_group = NEW_INTEGER(new_sz);
			for (int i = 0; i < state[4]; ++i)
				INTEGER(reads_per_group)[i] = orig_reads_per_group[i];
			SET_VECTOR_ELT(snap, READS_PER_GROUP_IDX, reads_per_group);
			state[4] = new_sz;
		}
		_trim_group_i(snap, state[1], state[2]);
		state[1] += 1;			/* next group */
		state[2] = 0;			/* reset offset */
		grp = NEW_RAW(GROUP_SIZE);
		SET_CURR_GRP(snap, grp);
	} 

	char *dest = (char *) RAW(grp);
	memcpy(dest + state[2], str, width);
	int i = state[0];
	INTEGER(VECTOR_ELT(snap, READS_PER_GROUP_IDX))[state[1]] += 1;
	INTEGER(VECTOR_ELT(snap, GROUP_IDX))[i] = state[1];
	INTEGER(VECTOR_ELT(snap, START_IDX))[i] = state[2] + 1;
	INTEGER(VECTOR_ELT(snap, WIDTH_IDX))[i] = width;
	state[0] += 1;
	state[2] += width;
}

SEXP
_to_XStringSet(SEXP seq, SEXP start, SEXP width, const char *baseclass,
			   const char *lkup)
{
	int seq_width = LENGTH(seq);

	char *str = (char *) RAW(seq);
	for (int i = 0; i < seq_width; ++i) {
		const char c = lkup[(int) str[i]];
		if (c == 0)
			Rf_error("invalid character '%c' in '%s'", str[i], baseclass);
		str[i] = c;
	}

	SEXP ptr = PROTECT(new_SharedVector("SharedRaw", seq));
	SEXP xstring = PROTECT(new_XVector(baseclass, ptr, 0, seq_width));
	SEXP irange = 
		PROTECT(new_IRanges("IRanges", start, width, R_NilValue));
	SEXP xstringset = new_XStringSet(NULL, xstring, irange);
	UNPROTECT(3);
	return xstringset;
}

const char *
_get_lookup(const char *baseclass)
{
	ENCODE_FUNC encode = encoder(baseclass);
	SEXP cls = PROTECT(eval(lang1(install(baseclass)), R_GlobalEnv));
	SEXP l = PROTECT(lang2(install("alphabet"), cls));
	SEXP alf = PROTECT(eval(l, R_GlobalEnv));

	char *lkup = (char *) R_alloc(256, sizeof(char));
	int i;
	if (alf == R_NilValue) {
		for (i = 0; i < 256; ++i) lkup[i] = (char) i;
	} else {
		for (i = 0; i < 256; ++i) lkup[i] = 0;
		for (i = 0; i < LENGTH(alf); ++i) {
			char c = CHAR(STRING_ELT(alf, i))[0];
			lkup[(int) c] = encode(c);
		}
		l = PROTECT(lang2(install("tolower"), alf));
		alf = PROTECT(eval(l, R_GlobalEnv));
		for (i = 0; i < LENGTH(alf); ++i) {
			char c = CHAR(STRING_ELT(alf, i))[0];
			lkup[(int) c] = encode(c);
		}
		UNPROTECT(2);
	}
	UNPROTECT(3);
	return lkup;
}

SEXP
_get_appender(const char *baseclass)
{
	char *class = (char *) R_alloc(strlen(baseclass) + 4, sizeof(char));
	sprintf(class, "%sSet", baseclass);
	SEXP l = PROTECT(lang3(install("selectMethod"), install("c"), 
						   mkString(class)));
	SEXP appender = eval(l, R_GlobalEnv);
	UNPROTECT(1);
	return appender;
}

SEXP
_XSnap_to_XStringSet(_XSnap snap)
{
	int *state = GET_STATE(snap);
	/* trim */
	_trim_group_i(snap, state[1], state[2]);

	int *reads_per_group = INTEGER(VECTOR_ELT(snap, READS_PER_GROUP_IDX));
	int *starts = INTEGER(VECTOR_ELT(snap, START_IDX));
	int *widths = INTEGER(VECTOR_ELT(snap, WIDTH_IDX));
	int tot_reads = 0;

	const char *baseclass = 
		CHAR(STRING_ELT(VECTOR_ELT(snap, BASENAME_IDX), 0));

	/* convert to xstringset; use list to manage protection */
	SEXP xstringset = PROTECT(NEW_LIST(state[1] + 1)); 
	const char *lkup = _get_lookup(baseclass);
	for (int i = 0; i <= state[1]; ++i) {
		SEXP start = PROTECT(NEW_INTEGER(reads_per_group[i]));
		SEXP width = PROTECT(NEW_INTEGER(reads_per_group[i]));
		int sz = reads_per_group[i] * sizeof(int);
		memcpy(INTEGER(start), starts + tot_reads, sz);
		memcpy(INTEGER(width), widths + tot_reads, sz);
		SEXP grp = GET_GRP_I(snap, i);
		SEXP xstringset1 = 
			_to_XStringSet(grp, start, width, baseclass, lkup);
		SET_VECTOR_ELT(xstringset, i, xstringset1);
		SET_GRP_I(snap, i, R_NilValue);
		UNPROTECT(2);
	}

	/* concatenate */
	SEXP appender = PROTECT(_get_appender(baseclass));
	int n = LENGTH(xstringset);
	while (n > 1) {
		SEXP res;
		int i;
		for (i = 0; i < n; i += 2) {
			if (i != n - 1) {
				SEXP l = lang3(appender, VECTOR_ELT(xstringset, i),
							   VECTOR_ELT(xstringset, i + 1));
				res = eval(l, R_GlobalEnv);
			} else {
				res = VECTOR_ELT(xstringset, i);
			}
			SET_VECTOR_ELT(xstringset, i, R_NilValue);
			SET_VECTOR_ELT(xstringset, i + 1, R_NilValue);
			SET_VECTOR_ELT(xstringset, i / 2, res);
		}
		n = i / 2;
	}
	
	UNPROTECT(2);
	return VECTOR_ELT(xstringset, 0);
}

void
_XSNAP_ELT(SEXP x, int elt)
{
	SEXP xstringset = 
		PROTECT(_XSnap_to_XStringSet(VECTOR_ELT(x, elt)));
	SET_VECTOR_ELT(x, elt, xstringset);
	UNPROTECT(1);
}
