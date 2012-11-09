#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include "ShortRead.h"
#include "IRanges_interface.h"

struct bufnode {
    int len;
    Rbyte *bytes;
    struct bufnode *next;
};

struct record {
    int length;
    const Rbyte *record;
};

struct records {
    int n, n_curr, n_tot, n_added;
    struct record *records;
};

SEXP _records_status(struct records *records)
{
    SEXP result = PROTECT(NEW_INTEGER(4));
    INTEGER(result)[0] = records->n;
    INTEGER(result)[1] = records->n_curr;
    INTEGER(result)[2] = records->n_added;
    INTEGER(result)[3] = records->n_tot;

    SEXP nms = PROTECT(NEW_CHARACTER(4));
    SET_STRING_ELT(nms, 0, mkChar("n"));
    SET_STRING_ELT(nms, 1, mkChar("current"));
    SET_STRING_ELT(nms, 2, mkChar("added"));
    SET_STRING_ELT(nms, 3, mkChar("total"));
    SET_NAMES(result, nms);

    UNPROTECT(2);
    return result;
}

/* fastq */

const Rbyte *_fastq_record_end(const Rbyte * buf, const Rbyte * bufend)
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

SEXP _fastq_as_XStringSet(struct records *fastq)
{
    SEXP widths = PROTECT(NEW_LIST(2));
    SET_VECTOR_ELT(widths, 0, NEW_INTEGER(fastq->n_curr));
    SET_VECTOR_ELT(widths, 1, NEW_INTEGER(fastq->n_curr));
    int *sread_w = INTEGER(VECTOR_ELT(widths, 0)),
        *id_w = INTEGER(VECTOR_ELT(widths, 1));

    /* geometry */
#ifdef SUPPORT_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < fastq->n_curr; ++i) {
        const Rbyte *buf = fastq->records[i].record;
        const Rbyte *start;

        start = ++buf;          /* id; skip '@' */
        while (!((*buf == '\n') || (*buf == '\r')))
            ++buf;
        id_w[i] = buf - start;
        while ((*buf == '\n') || (*buf == '\r'))
            ++buf;
        sread_w[i] = 0;         /* read */
        while (*buf != '+') {
            while (!((*buf == '\n') || (*buf == '\r'))) { /* strip '\n' */
                sread_w[i] += 1;
                ++buf;
            }
	    ++buf;
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

#ifdef SUPPORT_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < fastq->n_curr; ++i) {
        cachedCharSeq x;

        const Rbyte *buf = fastq->records[i].record,
            *bufend = buf + fastq->records[i].length,
            *start;
        char *curr;

        /* id */
        start = ++buf;          /* skip '@' */
        while (!((*buf == '\n') || (*buf == '\r')))
            ++buf;
        x = get_cachedXRawList_elt(&id, i);
        memcpy((char *) x.seq, start, (buf - start) * sizeof(Rbyte));

        /* read */
        while ((*buf == '\n') || (*buf == '\r'))
            ++buf;
        curr = (char *) get_cachedXRawList_elt(&sread, i).seq;
        while (*buf != '+') {
            while (!((*buf == '\n') || (*buf == '\r'))) /* strip '\n' */
                *curr++ = DNAencode(*buf++);
            buf++;
        }

        /* second id tag */
        while (!((*buf == '\n') || (*buf == '\r')))
            ++buf;

        /* quality */
        while ((*buf == '\n') || (*buf == '\r'))
            ++buf;              /* leading '\n' */
        start = buf;
        x = get_cachedXRawList_elt(&qual, i);
        curr = (char *) x.seq;
        while (buf != bufend && curr - x.seq != x.length) {
            if ((*buf != '\n') && (*buf != '\r'))
                *curr++ = *buf++;
            else
                buf++;
        }
    }
    SEXP nms = PROTECT(NEW_CHARACTER(3));
    SET_STRING_ELT(nms, 0, mkChar("sread"));
    SET_STRING_ELT(nms, 1, mkChar("quality"));
    SET_STRING_ELT(nms, 2, mkChar("id"));
    SET_NAMES(ans, nms);

    UNPROTECT(3);
    return ans;
}

/* Sampler */

struct sampler {
    struct records *sample;
    struct {
        struct record *records;
        int n, n_curr;
    } current;
    struct bufnode *bufnode;    /* tail end of binary stream */
};

struct sampler * _sampler_new(int n)
{
    struct sampler *sampler = Calloc(1, struct sampler);
    sampler->sample = Calloc(1, struct records);
    sampler->sample->records = Calloc(n, struct record);
    sampler->sample->n = n;
    sampler->current.records = Calloc(n, struct record);
    sampler->current.n = n;
    sampler->bufnode = Calloc(1, struct bufnode);
    return sampler;
}

void _sampler_reset(struct sampler *sampler)
{
    struct records *sample = sampler->sample;
    for (int i = 0; i < sample->n_curr; ++i)
        Free(sample->records[i].record);
    if (NULL != sampler->bufnode->bytes)
        Free(sampler->bufnode->bytes);
    sample->n_curr = sample->n_added = sample->n_tot = 0;
    sampler->current.n_curr = 0;
}

void _sampler_free(struct sampler *sampler)
{
    struct records *sample = sampler->sample;
    for (int i = 0; i < sample->n_curr; ++i)
        Free(sample->records[i].record);
    if (NULL != sampler->bufnode->bytes)
        Free(sampler->bufnode->bytes);
    Free(sampler->sample->records);
    Free(sampler->sample);
    Free(sampler->current.records);
    Free(sampler->bufnode);
    Free(sampler);
}

void _sampler_add1(struct records *sample, const Rbyte *record,
                   int len, int idx)
{
    /* add record to sample */
    if (sample->n_curr == sample->n)
        Free(sample->records[idx].record);

    sample->records[idx].length = len;
    Rbyte *intern_record = Calloc(len, Rbyte);
    memcpy(intern_record, record, len * sizeof(Rbyte));
    sample->records[idx].record = intern_record;
    sample->n_added += 1;
    sample->n_tot += 1;
}

int * _sampler_wout_replacement(int n, int k)
{
    /* sample k of n without replacement */
    int *idx = Calloc(n, int);
    for (int i = 0; i < n; ++i)
        idx[i] = i;
    for (int i = 0; i < k; ++i) {
        int j = (n - i) * unif_rand();
        int tmp = idx[i];
        idx[i] = idx[i + j];
        idx[i + j] = tmp;
    }
    return idx;
}

void _sampler_dosample(struct sampler *sampler)
{
    int n_curr = sampler->current.n_curr;
    int n_tot = n_curr + sampler->sample->n_tot;
    double n_choose = n_tot < sampler->sample->n ?
        n_tot : sampler->sample->n;
    int n_samp = rbinom(n_curr, n_choose / n_tot);

    if (0 != n_samp) {
        int sn_curr = sampler->sample->n_curr;
        int *keep = _sampler_wout_replacement(n_curr, n_samp);
        int *drop = _sampler_wout_replacement(sn_curr, n_samp);
        /* save selected reads */
        for (int i = 0; i < n_samp; ++i) {
            struct record *r = sampler->current.records + keep[i];
            _sampler_add1(sampler->sample, r->record, r->length, drop[i]);
        }

        Free(keep);
        Free(drop);
    }
    sampler->sample->n_tot = n_tot;
    sampler->current.n_curr = 0;
}

void _sampler_add(struct sampler *sampler, const Rbyte *record,
                  int len)
{
    struct records *sample = sampler->sample;
    if (sample->n_curr < sample->n) { /* sampling not yet needed */
        _sampler_add1(sample, record, len, sample->n_curr);
        sample->n_curr++;
    } else {                    /* sample */
        struct record *r =
            sampler->current.records + sampler->current.n_curr;
        r->record = record;
        r->length = len;
        if (sampler->current.n == ++sampler->current.n_curr)
            _sampler_dosample(sampler);
    }
}

void _sampler_scratch_set(struct sampler *sampler, const Rbyte *record,
                          int len)
{
    if (NULL != sampler->bufnode->bytes)
        Free(sampler->bufnode->bytes);
    if (NULL != record) {
        Rbyte *bytes = Calloc(len, Rbyte);
        memcpy(bytes, record, len * sizeof(Rbyte));
        sampler->bufnode->bytes = bytes;
    }
    sampler->bufnode->len = len;
}

/* R implementation -- FastqSampler */

#define SAMPLER(s) ((struct sampler *) R_ExternalPtrAddr(s))

void _sampler_finalize(SEXP s)
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
    R_RegisterCFinalizerEx(s, _sampler_finalize, TRUE);
    UNPROTECT(1);
    return s;
}

SEXP sampler_add(SEXP s, SEXP bin)
{
    /* create a buffer with both scratch and new data */
    struct sampler *sampler = SAMPLER(s);
    struct bufnode *scratch = sampler->bufnode;

    if (scratch->bytes) {
        int len = Rf_length(bin), buflen = scratch->len + len;
        Rbyte *buf = Calloc(buflen, Rbyte), *obuf = scratch->bytes;
        memcpy(buf, scratch->bytes, scratch->len * sizeof(Rbyte));
        Free(obuf);
        memcpy(buf + scratch->len, RAW(bin), len * sizeof(Rbyte));
        scratch->bytes = buf;
        scratch->len = buflen;
    } else {
        int buflen = Rf_length(bin);
        Rbyte *buf = Calloc(buflen, Rbyte);
        memcpy(buf, RAW(bin), buflen * sizeof(Rbyte));
        scratch->bytes = buf;
        scratch->len = buflen;
    }

    /* parse the buffer */
    const Rbyte *buf = scratch->bytes, *bufend = buf + scratch->len;
    GetRNGstate();
    while (buf < bufend) {
        while (buf < bufend && *buf == '\n')
            ++buf;
        const Rbyte *prev = buf;
        if (NULL == (buf = _fastq_record_end(buf, bufend))) {
            buf = prev;
            break;
        }
        _sampler_add(sampler, prev, buf - prev);
    }
    _sampler_dosample(sampler);
    PutRNGstate();

    if (bufend - buf) {
        int len = bufend - buf;
        Rbyte *tail = Calloc(len, Rbyte);
        memcpy(tail, buf, len * sizeof(Rbyte));
        Free(scratch->bytes);
        scratch->bytes = tail;
        scratch->len = len;
    } else {
        Free(scratch->bytes);
    }

    return s;
}

SEXP sampler_status(SEXP s)
{
    struct sampler *sampler = SAMPLER(s);
    return _records_status(sampler->sample);
}

SEXP sampler_as_XStringSet(SEXP s)
{
    struct sampler *sampler = SAMPLER(s);
    SEXP result = _fastq_as_XStringSet(sampler->sample);
    _sampler_scratch_set(sampler, NULL, 0);
    _sampler_reset(sampler);
    return result;
}

/* Streamer */

struct streamer {
    struct records *stream;
    struct bufnode *bufnode;
};

struct streamer * _streamer_new(int n)
{
    struct streamer *streamer = Calloc(1, struct streamer);
    streamer->stream = Calloc(1, struct records);
    streamer->stream->records = Calloc(n, struct record);
    streamer->stream->n = n;
    return streamer;
}

void _streamer_reset(struct streamer *streamer)
{
    streamer->stream->n_curr = 0;
    struct bufnode *bufnode = streamer->bufnode, *prev;
    if (NULL != bufnode) {
        bufnode = bufnode->next;
        while (NULL != bufnode) {
            prev = bufnode;
            bufnode = prev->next;
            Free(prev->bytes);
            Free(prev);
        }
        streamer->bufnode->next = NULL;
    }
}

void _streamer_free(struct streamer *streamer)
{
    struct bufnode *curr, *next = streamer->bufnode;
    while (next) {
        curr = next;
        next = curr->next;
        Free(curr->bytes);
        Free(curr);
    }
    Free(streamer->stream->records);
    Free(streamer->stream);
    Free(streamer);
}

void _streamer_add(struct records *stream, const Rbyte *record,
                   int len)
{
    stream->records[stream->n_curr].length = len;
    stream->records[stream->n_curr].record = record;
    stream->n_curr += 1;
    stream->n_added += 1;
}

#define STREAMER(s) ((struct streamer *) R_ExternalPtrAddr(s))

void _streamer_finalize(SEXP s)
{
    struct streamer *streamer = STREAMER(s);
    if (!streamer)
        return;
    _streamer_free(streamer);
    R_ClearExternalPtr(s);
}

SEXP streamer_new(SEXP n)
{
    struct streamer *streamer = _streamer_new(INTEGER(n)[0]);
    SEXP s = PROTECT(R_MakeExternalPtr(streamer, mkString("streamer"),
                                       R_NilValue));
    R_RegisterCFinalizerEx(s, _streamer_finalize, TRUE);
    UNPROTECT(1);
    return s;
}

SEXP streamer_add(SEXP s, SEXP bin, SEXP skipadd)
{
    struct streamer *streamer = STREAMER(s);
    int len = Rf_length(bin);
    int skip = INTEGER(skipadd)[0], add = INTEGER(skipadd)[1];

    /* start with tail of previous bin */
    struct bufnode *scratch = streamer->bufnode;
    if (NULL == scratch) {
        /* first record */
        scratch = streamer->bufnode = Calloc(1, struct bufnode);
    }
    if (NULL == scratch->bytes) {
        /* nothing 'extra' from previous bin */
        scratch->bytes = Calloc(len, Rbyte);
        scratch->len = len;
        memcpy(scratch->bytes, RAW(bin), len * sizeof(Rbyte));
    } else {
        /* scratch contains tail of prev. bin */
        int buflen = scratch->len;
        Rbyte *bytes = Calloc(buflen + len, Rbyte);
        memcpy(bytes, scratch->bytes, buflen * sizeof(Rbyte));
        memcpy(bytes + buflen, RAW(bin), len * sizeof(Rbyte));
        Free(scratch->bytes);
        scratch->bytes = bytes;
        scratch->len = buflen + len;
    }

    /* find record starts and lengths */
    const Rbyte *buf = scratch->bytes, *bufend = buf + scratch->len;
    struct records *stream = streamer->stream;
    while (add > stream->n_curr && buf < bufend) {
        while (buf < bufend && *buf == '\n')
            ++buf;
        const Rbyte *prev = buf;
        if (NULL == (buf = _fastq_record_end(buf, bufend))) {
            buf = prev;
            break;
        }

        stream->n_tot += 1;
        if (skip == 0)
            _streamer_add(stream, prev, buf - prev);
        else
            skip -= 1;
    }

    /* capture tail of bin */
    if (NULL != scratch->bytes) {
        struct bufnode *next = scratch;
        scratch = streamer->bufnode = Calloc(1, struct bufnode);
        scratch->next = next;
    }
    if (bufend - buf) {
        int len = bufend - buf;
        Rbyte *tail = Calloc(len, Rbyte);
        memcpy(tail, buf, len * sizeof(Rbyte));
        scratch->bytes = tail;
        scratch->len = len;
    }

    return s;
}

SEXP streamer_status(SEXP s)
{
    struct streamer *streamer = STREAMER(s);
    return _records_status(streamer->stream);
}

SEXP streamer_as_XStringSet(SEXP s)
{
    struct streamer *streamer = STREAMER(s);
    struct records *stream = streamer->stream;
    SEXP result = _fastq_as_XStringSet(stream);
    _streamer_reset(streamer);
    return result;
}
