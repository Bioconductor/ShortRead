/* 
 * An _XSNAP is a SEXP that contains sufficient information to create
 * (`snap') an XStringSet object from its content. It is allocated
 * once to an initial size, and grows as needed. Any `extra'
 * allocation is recovered.
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

static const int _BUFFERNODE_SIZE = 250000000;

/* _Buffer, _BufferNode: linked list of XString data chunks */

struct _Buffer {
    char *baseclass;
    int *offset, i_offset;
    struct _BufferNode *root, *curr;
};

struct _BufferNode {
    int n;                      /* number of entries */
    int buf_size;
    char *buf, *curr;
    /* linked list */
    struct _BufferNode *next;
};

/* _BufferNode implementation */

struct _BufferNode *_BufferNode_new()
{
    struct _BufferNode *node = Calloc(1, struct _BufferNode);
    node->curr = node->buf = Calloc(_BUFFERNODE_SIZE, char);
    node->n = 0;
    node->buf_size = _BUFFERNODE_SIZE;
    node->next = NULL;
    return node;
}

void _BufferNode_free(struct _BufferNode *node)
{
    Free(node->buf);
    Free(node);
}

void _BufferNode_encode(struct _BufferNode *node, const char *lkup)
{
    for (char *buf = node->buf; buf < node->curr; ++buf) {
        const char c = lkup[(int) *buf];
        if (c == 0)
            Rf_error("invalid character '%c'", c);
        *buf = c;
    }
}

int _BufferNode_append(struct _BufferNode *node, const char *s, int w)
{
    int offset = node->curr - node->buf;
    if (offset + w >= node->buf_size)
        return -1;
    memcpy(node->curr, s, w);
    node->curr += w;
    node->n += 1;
    return offset;
}

SEXP _BufferNode_snap(struct _BufferNode * node, const int *offset,
                      const char *baseclass)
{
    const int n_raw = node->curr - node->buf;
    SEXP seq = PROTECT(NEW_RAW(n_raw)),
        start = PROTECT(NEW_INTEGER(node->n)),
        width = PROTECT(NEW_INTEGER(node->n));
    memcpy(RAW(seq), node->buf, n_raw);
    for (int i = 0; i < node->n; ++i)
        INTEGER(start)[i] = offset[i] + 1;
    for (int i = 0; i < node->n - 1; ++i)
        INTEGER(width)[i] = offset[i + 1] - offset[i];
    if (node->n > 0)
        INTEGER(width)[node->n - 1] =
            node->curr - (node->buf + offset[node->n - 1]);
    SEXP xstringset = _to_XStringSet(seq, start, width, baseclass);
    UNPROTECT(3);
    return xstringset;
}

/* _Buffer implementation */

struct _Buffer *_Buffer_new(int n_offsets, const char *baseclass)
{
    struct _Buffer *buffer = Calloc(1, struct _Buffer);
    buffer->baseclass = Calloc(strlen(baseclass) + 1, char);
    buffer->offset = Calloc(n_offsets, int);
    buffer->i_offset = 0;
    strcpy(buffer->baseclass, baseclass);
    buffer->root = buffer->curr = _BufferNode_new();
    return buffer;
}

void _Buffer_free(struct _Buffer *buf)
{
    struct _BufferNode *curr = buf->root;
    while (curr != NULL) {
        struct _BufferNode *tmp = curr;
        curr = curr->next;
        _BufferNode_free(tmp);
    }
    Free(buf->offset);
    Free(buf->baseclass);
    Free(buf);
}

void _Buffer_append(struct _Buffer *buf, const char *s)
{
    int w = strlen(s);
    struct _BufferNode *curr = buf->curr;
    int i;
    if ((i = _BufferNode_append(curr, s, w)) < 0) {
        curr = buf->curr = curr->next = _BufferNode_new();
        i = _BufferNode_append(curr, s, w);
        if (i < 0)
            Rf_error("ShortRead internal: _BufferNode too small");
    }
    buf->offset[buf->i_offset++] = i;
}

void _Buffer_encode(struct _Buffer *buf)
{
    const char *lkup = _get_lookup(buf->baseclass);
    struct _BufferNode *curr;
    for (curr = buf->root; curr != NULL; curr = curr->next)
        _BufferNode_encode(curr, lkup);
}

SEXP _Buffer_snap(struct _Buffer *buf)
{
    int n_buf = 0, n_off = 0;
    struct _BufferNode *curr, *tmp;
    for (curr = buf->root; curr != NULL; curr = curr->next)
        ++n_buf;
    SEXP xstringsets = PROTECT(NEW_LIST(n_buf));
    curr = buf->root;
    for (int i = 0; i < n_buf; ++i) {
        SEXP xs = _BufferNode_snap(curr, buf->offset + n_off,
                                   buf->baseclass);
        SET_VECTOR_ELT(xstringsets, i, xs);
        n_off += curr->n;
        tmp = curr;
        curr = curr->next;
        _BufferNode_free(tmp);
    }
    buf->curr = buf->root = NULL;
    UNPROTECT(1);
    return xstringsets;
}

/* Wrap _Buffer in external pointer */
void _xsnap_finalizer(SEXP xsnap)
{
    struct _Buffer *buffer = R_ExternalPtrAddr(xsnap);
    if (!buffer)
        return;
    _Buffer_free(buffer);
    R_ClearExternalPtr(xsnap);
}

_XSnap _NEW_XSNAP(int n_elt, const char *baseclass)
{
    struct _Buffer *buffer = _Buffer_new(n_elt, baseclass);
    SEXP xsnap = PROTECT(R_MakeExternalPtr(buffer, PROTECT(mkString("XSnap")),
                                           R_NilValue));
    R_RegisterCFinalizerEx(xsnap, _xsnap_finalizer, TRUE);
    UNPROTECT(2);
    return xsnap;
}

void _APPEND_XSNAP(_XSnap snap, const char *str)
{
    _Buffer_append(R_ExternalPtrAddr(snap), str);
}

SEXP _to_XStringSet(SEXP seq, SEXP start, SEXP width, const char *baseclass)
{
    char classname[40];         /* longest string should be "DNAStringSet" */
    int res = snprintf(classname, sizeof(classname), "%sSet", baseclass);
    if (res < 0 || ((unsigned int) res) >= sizeof(classname))
        error("ShortRead internal error in _to_XStringSet(): "
              "'classname' buffer too small or other error");
    SEXP irange = PROTECT(new_IRanges("IRanges", start, width, R_NilValue));
    SEXP xstringset = new_XRawList_from_tag(classname, baseclass, seq, irange);
    UNPROTECT(1);
    return xstringset;
}

const char *_get_lookup(const char *baseclass)
{
    ENCODE_FUNC encode = encoder(baseclass);
    SEXP nmspc = PROTECT(_get_namespace("ShortRead"));
    SEXP lng1 = PROTECT(lang1(install(baseclass)));
    SEXP cls = PROTECT(eval(lng1, nmspc));
    SEXP lng2 = PROTECT(lang2(install("alphabet"), cls));
    SEXP alf = PROTECT(eval(lng2, nmspc));

    char *lkup = (char *) R_alloc(256, sizeof(char));
    int i;
    if (alf == R_NilValue) {
        for (i = 0; i < 256; ++i)
            lkup[i] = (char) i;
    } else {
        for (i = 0; i < 256; ++i)
            lkup[i] = 0;
        for (i = 0; i < LENGTH(alf); ++i) {
            char c = CHAR(STRING_ELT(alf, i))[0];
            lkup[(int) c] = encode(c);
        }
        lng2 = PROTECT(lang2(install("tolower"), alf));
        alf = PROTECT(eval(lng2, nmspc));
        for (i = 0; i < LENGTH(alf); ++i) {
            char c = CHAR(STRING_ELT(alf, i))[0];
            lkup[(int) c] = encode(c);
        }
        UNPROTECT(2);
    }
    UNPROTECT(5);
    return lkup;
}

SEXP _get_appender(const char *baseclass)
{
    char *class = (char *) R_alloc(strlen(baseclass) + 4, sizeof(char));
    sprintf(class, "%sSet", baseclass);
    SEXP cls = PROTECT(mkString(class));
    SEXP lng3 = PROTECT(lang3(install("selectMethod"), install("c"), cls));
    SEXP nmspc = PROTECT(_get_namespace("ShortRead"));
    SEXP appender = eval(lng3, nmspc);
    UNPROTECT(3);
    return appender;
}

SEXP _XSnap_to_XStringSet(_XSnap snap)
{
    struct _Buffer *buffer = (struct _Buffer *) R_ExternalPtrAddr(snap);
    _Buffer_encode(buffer);
    SEXP xstringset = PROTECT(_Buffer_snap(buffer));

    /* concatenate */
    SEXP appender = PROTECT(_get_appender(buffer->baseclass));
    SEXP nmspc = PROTECT(_get_namespace("ShortRead"));
    int n = LENGTH(xstringset);
    while (n > 1) {
        SEXP res;
        int i;
        for (i = 0; i < n; i += 2) {
            if (i != n - 1) {
                SEXP lng3 = PROTECT(lang3(appender, VECTOR_ELT(xstringset, i),
                                          VECTOR_ELT(xstringset, i + 1)));
                res = eval(lng3, nmspc);
                SET_VECTOR_ELT(xstringset, i + 1, R_NilValue);
                UNPROTECT(1);
            } else {
                res = VECTOR_ELT(xstringset, i);
            }
            SET_VECTOR_ELT(xstringset, i, R_NilValue);
            SET_VECTOR_ELT(xstringset, i / 2, res);
        }
        n = i / 2;
    }

    UNPROTECT(3);
    return VECTOR_ELT(xstringset, 0);
}

void _XSNAP_ELT(SEXP x, int elt)
{
    SEXP xstringset = PROTECT(_XSnap_to_XStringSet(VECTOR_ELT(x, elt)));
    SET_VECTOR_ELT(x, elt, xstringset);
    UNPROTECT(1);
}
