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

SEXP _to_XStringSet(SEXP seq, SEXP start, SEXP width, const char *baseclass);
const char *_get_lookup(const char *baseclass);

/* 50m 2x100 reads; 10G reads */
static const int _BUFFERNODE_SIZE = 100000000;
static const int _BUFFERNODE_OFFSET_SIZE = 5000000;

/* _Buffer, _BufferNode: linked list of XString data chunks */

struct _Buffer {
	char *baseclass;
	struct _BufferNode *root, *curr;
};

struct _BufferNode {
	int i_offset, n_offsets;
	int *offset;
	int buf_size;
	char *buf;
	/* linked list */
	struct _BufferNode *next;
};

/* _BufferNode implementation */

struct _BufferNode *
_BufferNode_new()
{
	struct _BufferNode *node = Calloc(1, struct _BufferNode);
	if (!node)
		Rf_error("ShortRead internal: failed to allocate _BufferNode");
	node->offset = Calloc(_BUFFERNODE_OFFSET_SIZE, int);
	if (!node->offset)
		Rf_error("ShortRead internal: failed to allocate _BufferNode offsets");
	node->offset[0] = 0;
	node->i_offset = 0;
	node->n_offsets = _BUFFERNODE_OFFSET_SIZE;
	node->buf = Calloc(_BUFFERNODE_SIZE, char);
	if (!node->buf)
		Rf_error("ShortRead internal: failed to allcoate _BufferNode buffer");
	node->buf_size = _BUFFERNODE_SIZE;
	node->next = 0;
	return node;
}

void
_BufferNode_free(struct _BufferNode *node)
{
	Free(node->buf);
	Free(node->offset);
	Free(node);
}

void
_BufferNode_encode(struct _BufferNode *node, const char *lkup)
{
	char *buf = node->buf;
	for (int i = 0; i < node->offset[node->i_offset]; ++i) {
		const char c = lkup[(int) buf[i]];
		if (c == 0)
			Rf_error("invalid character '%c'", c);
		buf[i] = c;
	}
}

void
_BufferNode_append(struct _BufferNode *node, const char *s, int w)
{
	memcpy(&node->buf[node->offset[node->i_offset]], s, w);
	node->i_offset += 1;
	node->offset[node->i_offset] = node->offset[node->i_offset-1] + w;
}

SEXP
_BufferNode_snap(struct _BufferNode *node, const char *baseclass)
{
	const int *offset = node->offset;
	SEXP seq = PROTECT(NEW_RAW(offset[node->i_offset])),
		start = PROTECT(NEW_INTEGER(node->i_offset)),
		width = PROTECT(NEW_INTEGER(node->i_offset));
	memcpy(RAW(seq), node->buf, offset[node->i_offset]);
	for (int i = 0; i < node->i_offset; ++i) {
		INTEGER(start)[i] = offset[i] + 1;
		INTEGER(width)[i] = offset[i + 1] - offset[i];
	}

	SEXP xstringset = _to_XStringSet(seq, start, width, baseclass);

	UNPROTECT(3);
	return xstringset;
}

/* _Buffer implementation */

struct _Buffer *
_Buffer_new(const char *baseclass)
{
	struct _Buffer *buffer = Calloc(1, struct _Buffer);
	if (!buffer)
		Rf_error("ShortRead internal: failed to allocate _Buffer");
	buffer->baseclass = Calloc(strlen(baseclass) + 1, char);
	strcpy(buffer->baseclass, baseclass);
	buffer->root = buffer->curr = _BufferNode_new();
	return buffer;
}

void
_Buffer_free(struct _Buffer *buf)
{
	struct _BufferNode *curr = buf->root;
	while (curr != NULL) {
		struct _BufferNode *tmp = curr;
		curr = curr->next;
		_BufferNode_free(tmp);
	}
	Free(buf->baseclass);
	Free(buf);
}

void
_Buffer_append(struct _Buffer *buf, const char *s)
{
	int w = strlen(s);
	struct _BufferNode *curr = buf->curr;
	if (curr->i_offset + 1 >= curr->n_offsets ||
		curr->offset[curr->i_offset] + w > curr->buf_size) {
		curr = buf->curr = curr->next = _BufferNode_new();
	}
	_BufferNode_append(curr, s, w);
}

void
_Buffer_encode(struct _Buffer *buf)
{
	const char *lkup = _get_lookup(buf->baseclass);
	struct _BufferNode *curr;
	for (curr = buf->root; curr != NULL; curr = curr->next)
		_BufferNode_encode(curr, lkup);
}

SEXP
_Buffer_snap(struct _Buffer *buf)
{
	int n_buf = 0;
	struct _BufferNode *curr;
	for (curr = buf->root; curr != NULL; curr = curr->next)
		++n_buf;
	SEXP xstringsets = PROTECT(NEW_LIST(n_buf));
	curr = buf->root;
	for (int i = 0; i < n_buf; ++i) {
		SET_VECTOR_ELT(xstringsets, i, 
					   _BufferNode_snap(curr, buf->baseclass));
		struct _BufferNode *tmp = curr;
		curr = curr->next;
		_BufferNode_free(tmp);
	}
	buf->curr = buf->root = NULL;
	UNPROTECT(1);
	return xstringsets;
}

/* Wrap _Buffer in external pointer */
void
_xsnap_finalizer(SEXP xsnap)
{
	struct _Buffer *buffer = R_ExternalPtrAddr(xsnap);
	if (!buffer) return;
	_Buffer_free(buffer);
	R_ClearExternalPtr(xsnap);
}

_XSnap
_NEW_XSNAP(int nelt, const char *baseclass)
{
	struct _Buffer *buffer = _Buffer_new(baseclass);
	SEXP xsnap = PROTECT(R_MakeExternalPtr(buffer, mkString("XSnap"),
										   R_NilValue));
	R_RegisterCFinalizerEx(xsnap, _xsnap_finalizer, TRUE);
	UNPROTECT(1);
	return xsnap;
}

void
_APPEND_XSNAP(_XSnap snap, const char *str)
{
	_Buffer_append(R_ExternalPtrAddr(snap), str);
}

SEXP
_to_XStringSet(SEXP seq, SEXP start, SEXP width, const char *baseclass)
{
	int seq_width = LENGTH(seq);

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
	struct _Buffer *buffer = (struct _Buffer *) R_ExternalPtrAddr(snap);
	_Buffer_encode(buffer);
	SEXP xstringset = PROTECT(_Buffer_snap(buffer));

	/* concatenate */
	SEXP appender = PROTECT(_get_appender(buffer->baseclass));
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
