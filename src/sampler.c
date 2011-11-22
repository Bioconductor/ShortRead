#include <stdlib.h>
#include "ShortRead.h"

#include "IRanges_interface.h"

/* fastq */

Rbyte *_fastq_record_end(Rbyte * buf, const Rbyte * bufend)
{
    int id = 1, nchr = 0;
    if (*buf++ != '@')
        Rf_error("record does not start with '@'");
    while (buf != bufend && *buf++ != '\n') ;	/* id 1 */
    while (buf != bufend && *buf != '+')	/* read */
        if (*buf++ != '\n')
            ++nchr;
    if (buf != bufend && *buf == '+')
        id -= 1;
    while (buf != bufend && *buf++ != '\n') ;	/* id 2 */
    while (buf != bufend && nchr)	/* qual */
        if (*buf++ != '\n')
            --nchr;
    if (0 != id || 0 != nchr)
        buf = NULL;
    if (buf && buf != bufend && *buf++ != '\n')
        Rf_error("internal: buf != <newline>");
    return buf;
}

/* sampler */

struct sampler_rec {
    int length;
    Rbyte *record;
};

struct sampler {
    int n, n_curr, n_tot, n_added;
    struct sampler_rec *records;
    struct sampler_rec *scratch;  /* tail end of binary stream */
};

struct sampler * _sampler_new(int n)
{
    struct sampler *sampler = Calloc(1, struct sampler);
    sampler->n = n;
    sampler->records = Calloc(n, struct sampler_rec);
    sampler->scratch = Calloc(1, struct sampler_rec);
    return sampler;
}

void _sampler_reset(struct sampler *sampler)
{
    for (int i = 0;i < sampler->n_curr; ++i)
        Free(sampler->records[i].record);
    sampler->n_tot = sampler->n_curr = sampler->n_added = 0;
    if (sampler->scratch->record)
        Free(sampler->scratch->record);
    sampler->scratch->record = NULL;
}

void _sampler_free(struct sampler *sampler)
{
    for (int i = 0; i < sampler->n_curr; ++i)
        Free(sampler->records[i].record);
    if (NULL != sampler->scratch->record)
        Free(sampler->scratch->record);
    Free(sampler->records);
    Free(sampler);
}

void _sampler_add(struct sampler *sampler, Rbyte *record, int len)
{
    int idx;
    sampler->n_tot += 1;

    if (sampler->n_curr < sampler->n) {
        idx = sampler->n_curr;
        sampler->n_curr++;
    } else {
        double r = (((double) rand()) / RAND_MAX);
        if (r >= ((double) sampler->n) / sampler->n_tot)
            return;
        idx = (((double) rand()) / RAND_MAX) * sampler->n;
        Free(sampler->records[idx].record);
    }
    sampler->n_added += 1;
    sampler->records[idx].length = len;
    sampler->records[idx].record = Calloc(len, Rbyte);
    memcpy(sampler->records[idx].record, record, len * sizeof(Rbyte));
}

void _sampler_scratch_set(struct sampler *sampler, Rbyte *record, int len)
{
    if (NULL != sampler->scratch->record)
        Free(sampler->scratch->record);
    if (NULL != record)
        sampler->scratch->record = Calloc(len, Rbyte);
    sampler->scratch->length = len;
    memcpy(sampler->scratch->record, record, len * sizeof(Rbyte));
}

struct sampler_rec * _sampler_scratch_get(struct sampler *sampler)
{
    return sampler->scratch;
}

SEXP _sampler_as_XStringSet(struct sampler *sampler)
{
    SEXP widths = PROTECT(NEW_LIST(2));
    SET_VECTOR_ELT(widths, 0, NEW_INTEGER(sampler->n_curr));
    SET_VECTOR_ELT(widths, 1, NEW_INTEGER(sampler->n_curr));
    int *sread_w = INTEGER(VECTOR_ELT(widths, 0)),
        *id_w = INTEGER(VECTOR_ELT(widths, 1));

    /* geometry */
#pragma omp parallel for
    for (int i = 0; i < sampler->n_curr; ++i) {
        const Rbyte *buf = sampler->records[i].record;
        const Rbyte *start;

        start = ++buf;          /* id; skip '@' */
        while (*buf != '\n')
            ++buf;
        id_w[i] = buf - start;
        while (*buf == '\n')
            ++buf;
        sread_w[i] = 0;         /* read */
        while (*buf != '+') {
            while (*buf++ != '\n') /* strip '\n' */
                sread_w[i] += 1;
        }
    }

    /* results */
    SEXP ans = PROTECT(NEW_LIST(3));
    SET_VECTOR_ELT(ans, 0, alloc_XRawList("DNAStringSet", "DNAString",
                                          VECTOR_ELT(widths, 0)));
    SET_VECTOR_ELT(ans, 1, alloc_XRawList("BStringSet", "BString",
                                          VECTOR_ELT(widths, 0)));
    SET_VECTOR_ELT(ans, 2, alloc_XRawList("BStringSet", "BString",
                                          VECTOR_ELT(widths, 1)));

    cachedXVectorList
        sread = cache_XVectorList(VECTOR_ELT(ans, 0)),
        qual = cache_XVectorList(VECTOR_ELT(ans, 1)),
        id = cache_XVectorList(VECTOR_ELT(ans, 2));

#pragma omp parallel for
    for (int i = 0; i < sampler->n_curr; ++i) {
        cachedCharSeq x;

        Rbyte *buf = sampler->records[i].record;
        Rbyte *bufend = buf + sampler->records[i].length;
        Rbyte *start, *curr;
        int nchr;

        /* id */
        start = ++buf;          /* skip '@' */
        while (*buf != '\n')
            ++buf;
        x = get_cachedXRawList_elt(&id, i);
        memcpy((char *) x.seq, start, (buf - start) * sizeof(Rbyte));

        /* read */
        while (*buf == '\n')
            ++buf;
        start = curr = buf;
        while (*buf != '+') {
            while (*buf != '\n') /* strip '\n' */
                *curr++ = DNAencode(*buf++);
            buf++;
        }
        x = get_cachedXRawList_elt(&sread, i);
        memcpy((char *) x.seq, start, (curr - start) * sizeof(Rbyte));
        nchr = curr - start;

        /* second id tag */
        while (*buf != '\n')
            ++buf;

        /* quality */
        while (*buf == '\n')
            ++buf;              /* leading '\n' */
        start = curr = buf;
        while (buf != bufend && curr - start != nchr) {
            if (*buf != '\n')
                *curr++ = *buf++;
            else
                buf++;
        }
        x = get_cachedXRawList_elt(&qual, i);
        memcpy((char *) x.seq, start, (curr - start) * sizeof(Rbyte));
    }

    _sampler_reset(sampler);

    SEXP nms = PROTECT(NEW_CHARACTER(3));
    SET_STRING_ELT(nms, 0, mkChar("sread"));
    SET_STRING_ELT(nms, 1, mkChar("quality"));
    SET_STRING_ELT(nms, 2, mkChar("id"));
    SET_NAMES(ans, nms);

    UNPROTECT(3);
    return ans;
}

/* R implementation */

#define SAMPLER(s) ((struct sampler *) R_ExternalPtrAddr(s))
#define SCRATCH(s) (_sampler_scratch_get(SAMPLER(s)))

void _sampler_finalizer(SEXP s)
{
    struct sampler *sampler = SAMPLER(s);
    if (!sampler)
        return;
    _sampler_free(sampler);
    R_ClearExternalPtr(s);
}

SEXP sampler_new(SEXP n)
{
    struct sampler *sampler = _sampler_new(INTEGER(n)[0]);
    SEXP s = PROTECT(R_MakeExternalPtr(sampler, mkString("sampler"),
                                       R_NilValue));
    R_RegisterCFinalizerEx(s, _sampler_finalizer, TRUE);
    UNPROTECT(1);
    return s;
}

SEXP sampler_add(SEXP s, SEXP bin)
{
    /* create a buffer with both scratch and new data */
    struct sampler_rec *scratch = SCRATCH(s);
    int buflen = scratch->length + Rf_length(bin);
    Rbyte *buf = (Rbyte *) R_alloc(sizeof(Rbyte), buflen);
    memcpy(buf, scratch->record, scratch->length * sizeof(Rbyte));
    memcpy(buf + scratch->length, RAW(bin),
           Rf_length(bin) * sizeof(Rbyte));

    /* parse the buffer */
    struct sampler *sampler = SAMPLER(s);
    Rbyte *prev = buf, *bufend = buf + buflen;
    _sampler_scratch_set(sampler, NULL, 0);
    while (buf && buf < bufend) {
        while (buf < bufend && *buf == '\n')
            ++buf;
        prev = buf;
        if (NULL == (buf = _fastq_record_end(buf, bufend))) {
            _sampler_scratch_set(sampler, prev, bufend - prev);
            break;
        }
        _sampler_add(sampler, prev, buf - prev);
    }

    return s;
}

SEXP sampler_summary(SEXP s)
{
    struct sampler *sampler = SAMPLER(s);
    SEXP result = PROTECT(NEW_INTEGER(4));
    INTEGER(result)[0] = sampler->n;
    INTEGER(result)[1] = sampler->n_curr;
    INTEGER(result)[2] = sampler->n_added;
    INTEGER(result)[3] = sampler->n_tot;

    SEXP nms = PROTECT(NEW_CHARACTER(4));
    SET_STRING_ELT(nms, 0, mkChar("n"));
    SET_STRING_ELT(nms, 1, mkChar("current"));
    SET_STRING_ELT(nms, 2, mkChar("added"));
    SET_STRING_ELT(nms, 3, mkChar("total"));
    SET_NAMES(result, nms);

    UNPROTECT(2);
    return result;
}

SEXP sampler_as_XStringSet(SEXP s)
{
    return _sampler_as_XStringSet(SAMPLER(s));
}

/*  */

SEXP sampler_rec_counter(SEXP buffer)
{
    /* full and partial records */
    if (0 == LENGTH(buffer))
        return ScalarInteger(0);

    Rbyte *buf = RAW(buffer);
    const Rbyte *bufend = buf + LENGTH(buffer);
    int n = 0;

    while (buf && buf < bufend) {
        while (buf != bufend && *buf == '\n')
            ++buf;
        if (buf != bufend) {
            buf = _fastq_record_end(buf, bufend);
            ++n;
        }
    }
    return ScalarInteger(n);
}

SEXP sampler_rec_parser(SEXP buffer, SEXP rec_n)
{
    if (0 == INTEGER(rec_n)[0])
        return NEW_LIST(0);

    const int n = INTEGER(rec_n)[0];
    Rbyte *buf = RAW(buffer);
    const Rbyte *bufend = buf + LENGTH(buffer), *prev;
    SEXP lst, elt;
    int s_idx;

    lst = PROTECT(NEW_LIST(n));

    prev = buf = RAW(buffer);
    s_idx = 0;
    while (buf && buf < bufend) {
        while (buf < bufend && *buf == '\n')
            ++buf;
        prev = buf;
        if (NULL == (buf = _fastq_record_end(buf, bufend)))
            break;
        elt = NEW_RAW(buf - prev);
        SET_VECTOR_ELT(lst, s_idx++, elt);
        memcpy((Rbyte *) RAW(elt), prev, (buf - prev) * sizeof(Rbyte));
    }
    if (n > s_idx) {
        elt = NEW_RAW(bufend - prev);
        SET_VECTOR_ELT(lst, s_idx++, elt);
        memcpy((Rbyte *) RAW(elt), prev, (bufend - prev) * sizeof(Rbyte));
    }
    if (n != s_idx)
        Rf_error("internal: sample index (%d) != 'rec_n' (%d)",
                 s_idx, n);
    UNPROTECT(1);
    return lst;
}

SEXP _as_xstringset(SEXP width, SEXP tag, const char *classname,
                    const char *element_type)
{
    int nrec = LENGTH(width);
    SEXP start = PROTECT(NEW_INTEGER(nrec));
    int i, *s = INTEGER(start);
    const int *w = INTEGER(width);
    s[0] = 1;
    for (i = 1; i < nrec; ++i)
        s[i] = s[i - 1] + w[i - 1];
    SEXP rng = PROTECT(new_IRanges("IRanges", start, width, R_NilValue));
    SEXP xstringset =
        PROTECT(new_XRawList_from_tag(classname, element_type, tag, rng));
    UNPROTECT(3);
    return xstringset;
}

SEXP sampler_as_fastq(SEXP records)
{
    SEXP ans = PROTECT(NEW_LIST(3));
    SEXP lengths = PROTECT(NEW_LIST(2));
    int nrec = LENGTH(records), totlen = 0, i;

    SET_VECTOR_ELT(lengths, 0, NEW_INTEGER(nrec));	/* sread / quality */
    SET_VECTOR_ELT(lengths, 1, NEW_INTEGER(nrec));	/* id */

    for (i = 0; i < nrec; ++i)
        totlen += LENGTH(VECTOR_ELT(records, i));
    SET_VECTOR_ELT(ans, 0, NEW_RAW(totlen / 2));	/* sread */
    SET_VECTOR_ELT(ans, 1, NEW_RAW(totlen / 2));	/* quality */
    SET_VECTOR_ELT(ans, 2, NEW_RAW(totlen / 2));	/* id */

    Rbyte *sread_offset = RAW(VECTOR_ELT(ans, 0)),
        *qual_offset = RAW(VECTOR_ELT(ans, 1)),
        *id_offset = RAW(VECTOR_ELT(ans, 2));
    int *sread_len = INTEGER(VECTOR_ELT(lengths, 0)),
        *id_len = INTEGER(VECTOR_ELT(lengths, 1));

    for (i = 0; i < nrec; ++i) {
        SEXP record = VECTOR_ELT(records, i);
        Rbyte *buf = RAW(record);
        Rbyte *bufend = buf + LENGTH(record);
        Rbyte *start, *curr;
        int nchr;

        /* id */
        start = ++buf;          /* skip '@' */
        while (*buf != '\n')
            ++buf;
        memcpy(id_offset, start, (buf - start) * sizeof(Rbyte));
        *id_len++ = buf - start;
        id_offset += buf - start;

        /* read */
        while (*buf == '\n')
            ++buf;
        start = curr = buf;
        while (*buf != '+') {
            /* strip '\n' */
            while (*buf != '\n')
                *curr++ = *buf++;
            buf++;
        }
        memcpy(sread_offset, start, (curr - start) * sizeof(Rbyte));
        *sread_len++ = curr - start;
        sread_offset += curr - start;
        nchr = curr - start;

        /* second id tag */
        while (*buf != '\n')
            ++buf;
        /* quality */
        while (*buf == '\n')
            ++buf;              /* leading '\n' */
        start = curr = buf;
        while (buf != bufend && curr - start != nchr) {
            if (*buf != '\n')
                *curr++ = *buf++;
            else
                buf++;
        }
        memcpy(qual_offset, start, (curr - start) * sizeof(Rbyte));
        qual_offset += curr - start;
    }

    Rbyte *dna = RAW(VECTOR_ELT(ans, 0));
    while (dna < sread_offset) {
        Rbyte tmp = DNAencode(*dna);
        *dna++ = tmp;
    }

    SETLENGTH(VECTOR_ELT(ans, 0), sread_offset - RAW(VECTOR_ELT(ans, 0)));
    SETLENGTH(VECTOR_ELT(ans, 1), qual_offset - RAW(VECTOR_ELT(ans, 1)));
    SETLENGTH(VECTOR_ELT(ans, 2), id_offset - RAW(VECTOR_ELT(ans, 2)));

    SEXP xsset;
    xsset = _as_xstringset(VECTOR_ELT(lengths, 0), VECTOR_ELT(ans, 0),
                           "DNAStringSet", "DNAString");
    SET_VECTOR_ELT(ans, 0, xsset);;
    xsset = _as_xstringset(VECTOR_ELT(lengths, 0), VECTOR_ELT(ans, 1),
                           "BStringSet", "BString");
    SET_VECTOR_ELT(ans, 1, xsset);;
    xsset = _as_xstringset(VECTOR_ELT(lengths, 1), VECTOR_ELT(ans, 2),
                           "BStringSet", "BString");
    SET_VECTOR_ELT(ans, 2, xsset);

    SEXP nms = PROTECT(NEW_CHARACTER(3));
    SET_STRING_ELT(nms, 0, mkChar("sread"));
    SET_STRING_ELT(nms, 1, mkChar("quality"));
    SET_STRING_ELT(nms, 2, mkChar("id"));
    SET_NAMES(ans, nms);

    UNPROTECT(3);
    return ans;
}
